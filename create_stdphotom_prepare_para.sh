#!/bin/bash -u

# create catalogues suitable for absolute photometric calibration with
# 'create_abs_photo_info.sh'. The script takes astrometric information
# from a scamp processing and calculates 'correct' sky coordinates
# from objects catalogues (extracted from standard star
# fields). Afterwards the standard star object catalogues are merged
# with a reference standard star sources (e.g. from Landolt, Stetson,
# Sloan ....)

# SCRIPT HISTORY:
#
# 23.11.2012:
# The merging of extarcted source lists and the standard star catalogue
# needed an update because names of magnitudes are different for
# different filter systems. They are now read dynamically from the
# standardstar catalogue and a config file for 'make_ssc' is created on
# the fly.
#
# 07.01.2013:
# I made the script more robust to non-existing files.

#$1: main directory
#$2: science dir. (the catalogues are assumed to be in $1/$2/cat)
#$3: image extension
#$4: astrometry standardstar catalogue used for the scamp calibration
#    (the script needs the scamp headers which are assumed to be in 
#    $1/$2/headers_scamp_$4)
#$5: photometric standard star catalogue with full path
#$6: chips to work on

. ./${INSTRUMENT:?}.ini
. ./bash_functions.include
theli_start "$*" "${!#}"

# check number of command line arguments:
if [ $# -ne 6 ]; then
  theli_error "Number of command line argument not correct!" "${!#}"
  exit 1;
fi

# give meaningful names to command line arguments:
MD=$1
SD=$2
EXTENSION=$3
STANDARDCAT=$4
PHOTCAT=$5

# create configuration file 'make_ssc' (merging of extracted sources
# with standard star catalogue below) on the fly:
grep "COL_INPUT"  ${DATACONF}/stdphotom_prepare_make_ssc.conf |\
  ${P_GAWK} '{print $3}' > confraw.txt_$$

echo "SeqNr" >> confraw.txt_$$

# get keys that need to be added to stdphotom_prepare_make_ssc.conf
# from the standardstar catalogue. For simplicity we consider all keys
# present in the standardstar catalogue but not in our source lists:
ldacdesc -i ${PHOTCAT} -t STDTAB |\
   grep "Key name:" | ${P_GAWK} -F. '{print $NF}' > photcat.txt_$$

cat confraw.txt_$$  photcat.txt_$$ | sort | uniq -c |\
   ${P_GAWK} '($1 == 1) {print $2}' > uniq.txt_$$

cat  photcat.txt_$$ uniq.txt_$$ | sort | uniq -c |\
   ${P_GAWK} '($1 == 2) {print $2}' > add.txt_$$

cat ${DATACONF}/stdphotom_prepare_make_ssc.conf > \
  ${TEMPDIR}/make_ssc.conf_$$

# keys that are added from the second (= 1) catalogue
while read COL
do
  {
    echo "#"
    echo "COL_NAME  = ${COL}"
    echo "COL_INPUT = ${COL}"
    echo "COL_MERGE = AVE_REG"
    echo "COL_CHAN  = 1"
  } >>  ${TEMPDIR}/make_ssc.conf_$$
done < add.txt_$$

# keys that are added from the first (= 0) catalogue
echo "Xpos" > add2.txt_$$
echo "Ypos" >> add2.txt_$$
echo "Xpos_global" >> add2.txt_$$
echo "Ypos_global" >> add2.txt_$$

while read COL
do
  {
    echo "#"
    echo "COL_NAME  = ${COL}"
    echo "COL_INPUT = ${COL}"
    echo "COL_MERGE = AVE_REG"
    echo "COL_CHAN  = 0"
  } >>  ${TEMPDIR}/make_ssc.conf_$$
done < add2.txt_$$


# do the real work now:
for CHIP in ${!#}
do
  CATS=`find ${MD}/${SD}/cat -name \*_${CHIP}${EXTENSION}.cat`

  if [ "${CATS}" != "" ]; then
    for CAT in ${CATS}
    do
      BASE=`basename ${CAT} .cat`
      HEADBASE=`basename ${CAT} ${EXTENSION}.cat`
      ${P_LDACCONV} -i ${CAT} -o ${MD}/${SD}/cat/${BASE}_ldac.cat \
                    -b ${CHIP} -c ${INSTRUMENT} -f dummy 
  
      ${P_LDACTOASC} -b -i ${MD}/${SD}/cat/${BASE}_ldac.cat \
           -t OBJECTS -k ALPHA_J2000 DELTA_J2000 > ${TEMPDIR}/${BASE}_ldac.asc
  
      # get 'correct' astrometric scamp information to the catalogues:    
      ${S_APLASTROMSCAMP} -i ${MD}/${SD}/cat/${BASE}_ldac.cat \
           -o ${MD}/${SD}/cat/${BASE}_ldac_corr.cat \
           -s ${MD}/${SD}/headers_scamp_${STANDARDCAT}/${HEADBASE}.head \
           -t OBJECTS -p Xpos Ypos -r Ra Dec

      # get global xy coodinates with respect to the middle of the camera
      ${P_SKY2XY} -x 0.0 0.0 ${MD}/${SD}/${BASE}.fits @${TEMPDIR}/${BASE}_ldac.asc | ${P_GAWK} '{print $5, $6}' \
           > ${TEMPDIR}/${BASE}_ldac_global.asc
      ${P_ASCTOLDAC} -i ${TEMPDIR}/${BASE}_ldac_global.asc \
           -o ${TEMPDIR}/${BASE}_ldac_global.cat \
           -t OBJECTS -c ${DATACONF}/stdphotom_prepare_global_coordinates.conf
      ${P_LDACJOINKEY} -i ${MD}/${SD}/cat/${BASE}_ldac_corr.cat \
           -p ${TEMPDIR}/${BASE}_ldac_global.cat \
           -o ${MD}/${SD}/cat/${BASE}_ldac_corr_global_coordinates.cat \
           -t OBJECTS -k Xpos_global Ypos_global
      rm ${TEMPDIR}/${BASE}_ldac.asc ${TEMPDIR}/${BASE}_ldac_global.asc ${TEMPDIR}/${BASE}_ldac_global.cat
  
      ${P_MAKEJOIN} -i /${MD}/${SD}/cat/${BASE}_ldac_corr_global_coordinates.cat \
                    -o ${TEMPDIR}/tmp.cat0_$$ \
  		  -c ${DATACONF}/stdphotom_prepare_make_join.conf
      # calculate object magnitudes with a magzeropoint of 0 and
      # an exposure time normalisation of 1s (1.08574 is 2.5 / log(10)):
      ${P_LDACCALC} -i ${TEMPDIR}/tmp.cat0_$$ \
                    -o ${TEMPDIR}/tmp.cat00_$$ \
                    -t OBJECTS -c "(-1.08574*log(FLUX_AUTO/EXPTIME));" \
                    -n MAG_AUTO_corr "exposure time corr. MAG_AUTO" -k FLOAT
      ${P_LDACRENTAB} -i ${TEMPDIR}/tmp.cat00_$$\
                      -o ${TEMPDIR}/tmp.cat1_$$ -t OBJECTS STDTAB
      ${P_LDACRENKEY} -i ${TEMPDIR}/tmp.cat1_$$ -o ${TEMPDIR}/tmp.cat2_$$ \
                      -t STDTAB -k A_WORLD A_WCS B_WORLD B_WCS \
                                   THETA_J2000 THETAWCS
      ${P_ASSOCIATE} -i ${TEMPDIR}/tmp.cat2_$$ ${PHOTCAT}\
  		   -o ${TEMPDIR}/tmp.cat3_$$ ${TEMPDIR}/tmp.cat4_$$ \
                     -t STDTAB -c ${DATACONF}/stdphotom_prepare_associate.conf
  
      ${P_LDACFILTER} -i ${TEMPDIR}/tmp.cat3_$$ -o ${TEMPDIR}/tmp.cat5_$$ \
                      -c "(Pair_1>0);" -t STDTAB
  
      if [ "$?" -eq "0" ]; then
        ${P_LDACFILTER} -i ${TEMPDIR}/tmp.cat4_$$ \
                        -o ${TEMPDIR}/tmp.cat6_$$ -c "(Pair_0>0);" -t STDTAB
        ${P_ASSOCIATE} -i ${TEMPDIR}/tmp.cat5_$$ ${TEMPDIR}/tmp.cat6_$$ \
                       -o ${TEMPDIR}/tmp.cat7_$$ ${TEMPDIR}/tmp.cat8_$$ \
        	       -t STDTAB -c ${DATACONF}/stdphotom_prepare_associate.conf
        ${P_MAKESSC} -i ${TEMPDIR}/tmp.cat7_$$ ${TEMPDIR}/tmp.cat8_$$ \
           	   -o /${MD}/${SD}/cat/${BASE}_merg.cat\
          	   -t STDTAB \
                     -c ${TEMPDIR}/make_ssc.conf_$$
      fi
    done
    ${P_LDACPASTE} -i /${MD}/${SD}/cat/*_${CHIP}${EXTENSION}_merg.cat \
                   -o /${MD}/${SD}/cat/chip_${CHIP}_merg.cat -t PSSC
  else # if [ "${CATS}" != "" ]
    theli_warn "No catalogues for Chip ${CHIP} available!" "${!#}"
  fi
done

# clean up and bye
rm  ${TEMPDIR}/*_$$

theli_end "${!#}"
exit 0;
