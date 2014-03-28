#!/bin/bash -u

# The script merges astrometrically calbrated standardstar observations
# with a standardstar catalogue. It represents the second of a necessary
# to-stage process. This script cannot be used without having successfully
# run 'create_stdphotom_prepare_para.sh' before.

# SCRIPT HISTORY:
# ===============
#
# 05.06.2013:
# script started

#$1: main directory
#$2: science dir. (the catalogues are assumed to be in $1/$2/cat)
#$3: image extension
#$4: photometric standard star catalogue with full path

. ./${INSTRUMENT:?}.ini
. ./bash_functions.include
theli_start "$*"

# check number of command line arguments:
if [ $# -ne 4 ]; then
  theli_error "Number of command line argument not correct!"
  exit 1;
fi

# give meaningful names to command line arguments:
MD=$1
SD=$2
EXTENSION=$3
PHOTCAT=$4

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
echo "BADCCD" >> add2.txt_$$

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

# get a list of all the exposures:
${P_FIND} ${MD}/${SD}/cat -name \*${EXTENSION}_photprep.cat |\
  ${P_GAWK} -F/ '{print $NF}' | sed -e 's/'${EXTENSION}'_photprep\.cat//' | \
  ${P_GAWK} -F_ '{$NF = ""; $1 = $1; print $0}' | sort | uniq >  \
    ${TEMPDIR}/exposures.txt_$$

# do we need to do something at all?
if [ -s ${TEMPDIR}/exposures.txt_$$ ]; then

  # distribute task on different processors/cores:
  NIMA=`wc -l ${TEMPDIR}/exposures.txt_$$ | ${P_GAWK} '{print $1}'`
  
  if [ ${NPARA} -gt ${NIMA} ]; then
    NPARA=${NIMA}
  fi
  
  PROC=0
  while [ ${PROC} -lt ${NPARA} ]
  do
    NFIELDS[${PROC}]=0
    FIELDS[${PROC}]=""
    PROC=$(( ${PROC} + 1 ))
  done
  
  i=1
  while read FIELD
  do
    PROC=$(( $i % ${NPARA} ))
    FIELDS[${PROC}]="${FIELDS[${PROC}]} ${FIELD}"
    NFIELDS[${PROC}]=$(( ${NFIELDS[${PROC}]} + 1 ))
    i=$(( $i + 1 ))
  done < ${TEMPDIR}/exposures.txt_$$
  
  PROC=0
  while [ ${PROC} -lt ${NPARA} ]
  do
    j=$(( ${PROC} + 1 )) 
    echo -e "Starting Job $j. It has ${NFIELDS[${PROC}]} files to process!\n" 
    {
      k=1
      for FIELD in ${FIELDS[${PROC}]}
      do     
        echo "Working on exposure ${FIELD}"
        ${P_LDACPASTE} -i /${MD}/${SD}/cat/${FIELD}*${EXTENSION}_photprep.cat \
                       -o ${TEMPDIR}/tmp.cat0_${j}_$$ -t STDTAB

        # note that ldacpaste does not update SeqNr. So we have to do this
        # manually:
        ${P_LDACDELKEY} -i ${TEMPDIR}/tmp.cat0_${j}_$$ \
                        -o ${TEMPDIR}/tmp.cat1_${j}_$$ \
                        -t STDTAB -k SeqNr

        ${P_LDACADDKEY} -i ${TEMPDIR}/tmp.cat1_${j}_$$ \
                        -o /${MD}/${SD}/cat/${FIELD}_all_photprep.cat \
                        -t STDTAB -k SeqNr 1 COUNT "running object number"

        # with the correct SeqNrs we can finally associate:
        ${P_ASSOCIATE} -i /${MD}/${SD}/cat/${FIELD}_all_photprep.cat \
                          ${PHOTCAT}\
		       -o ${TEMPDIR}/tmp.cat3_${j}_$$ \
                          ${TEMPDIR}/tmp.cat4_${j}_$$ \
                       -t STDTAB \
                       -c ${DATACONF}/stdphotom_prepare_associate.conf

        ${P_MAKESSC} -i ${TEMPDIR}/tmp.cat3_${j}_$$ \
                        ${TEMPDIR}/tmp.cat4_${j}_$$ \
       	             -o ${TEMPDIR}/tmp.cat5_${j}_$$ -t STDTAB \
                     -c ${TEMPDIR}/make_ssc.conf_$$

        ${P_LDACFILTER} -i ${TEMPDIR}/tmp.cat5_${j}_$$ \
                        -o /${MD}/${SD}/cat/${FIELD}_all_photprep_merg.cat\
                        -c "((N_00=1)AND(N_01=1));" -t PSSC
  
        k=$(( $k + 1 ))
      done
    } &
    PROC=$(( ${PROC} + 1 ))
  done
  
  # only continue once all processes have finished
  wait;
else
  theli_error "No images to process!"
  cleanTmpFiles
  exit 1;
fi

# clean up and bye:
rm  ${TEMPDIR}/*_$$

theli_end
exit 0;
