#!/bin/bash -u

# First step to create catalogues suitable for absolute photometric
# calibration with 'create_abs_photo_info.sh'. The script takes
# astrometric information from a scamp processing and calculates
# 'correct' sky coordinates from objects catalogues (extracted from
# standard star fields). Afterwards the standard star object
# catalogues are merged with a reference standard star sources.
# (e.g. from Landolt, Stetson, Sloan ....)
# This last step is to be performed with the script
# 'create_stadphotom_merge.sh'

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
#
# 04.06.2013:
# I completely rewrote this script by splitting it. The task of estimating
# correct astrometry for standard star observations and merging them
# with a standard star catalogue is now donw in two separate scripts.
# The merging is no longer performed chip-by-chip but on an exposure level.
# This dramatically reduces the number of necessary associations. The
# associations became a very long process when we used also KiDS science
# data (overlapping with SDSS) in the calibration process.

#$1: main directory
#$2: science dir. (the catalogues are assumed to be in $1/$2/cat)
#$3: image extension
#$4: astrometry standardstar catalogue used for the scamp calibration
#    (the script needs the scamp headers which are assumed to be in 
#    $1/$2/headers_scamp_$4)
#$5: chips to work on

. ./${INSTRUMENT:?}.ini
. ./bash_functions.include
theli_start "$*" "${!#}"

# check number of command line arguments:
if [ $# -ne 5 ]; then
  theli_error "Number of command line argument not correct!" "${!#}"
  exit 1;
fi

# give meaningful names to command line arguments:
MD=$1
SD=$2
EXTENSION=$3
STANDARDCAT=$4

# start script task:
for CHIP in ${!#}
do
  ALLCATS=`find ${MD}/${SD}/cat -name \*_${CHIP}${EXTENSION}.cat`

  if [ "${ALLCATS}" != "" ]; then
    FIRSTCAT=`echo ${ALLCATS} | ${P_GAWK} '{print $1}'`
    FIRSTFILE=`basename ${FIRSTCAT} .cat`
    cp ${MD}/${SD}/${FIRSTFILE}.fits ${TEMPDIR}/dummy.fits_$$

    CATS=`echo`
    NUMBAD=0
    NUMOK=0
    for CAT in ${ALLCATS}
    do
      BASE=`basename ${CAT} .cat`
      BADCCD=`${P_DFITS} /${MD}/${SD}/${BASE}.fits | fitsort BADCCD | grep ${BASE} | ${P_GAWK} '{print $2}'`
	if [ "${BADCCD}" == "0" ]; then
	  CATS=`echo ${CATS} ${CAT}`
	  NUMOK=$(( NUMOK + 1 ))
	else
	  NUMBAD=$(( NUMBAD + 1 ))
	fi
    done

    if [ "${NUMBAD}" != "0" ]; then
      theli_warning "${NUMBAD}/$(( NUMOK + NUMBAD )) files of chip ${CHIP} have a BADCCD flag."
    fi

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
      CRVAL1=`grep CRVAL1 ${MD}/${SD}/headers_scamp_${STANDARDCAT}/${HEADBASE}.head | ${P_GAWK} '{print $3}'`
      value ${CRVAL1}
      writekey ${TEMPDIR}/dummy.fits_$$ CRVAL1 "${VALUE} / WCS Ref. value (RA in decimal degrees)" REPLACE

      CRVAL2=`grep CRVAL2 ${MD}/${SD}/headers_scamp_${STANDARDCAT}/${HEADBASE}.head | ${P_GAWK} '{print $3}'`
      value ${CRVAL2}
      writekey ${TEMPDIR}/dummy.fits_$$ CRVAL2 "${VALUE} / WCS Ref. value (DEC in decimal degrees)" REPLACE

      CRPIX1=`grep CRPIX1 ${MD}/${SD}/headers_scamp_${STANDARDCAT}/${HEADBASE}.head | ${P_GAWK} '{print $3}'`
      value ${CRPIX1}
      writekey ${TEMPDIR}/dummy.fits_$$ CRPIX1 "${VALUE} / WCS Coordinate reference pixel" REPLACE

      CRPIX2=`grep CRPIX2 ${MD}/${SD}/headers_scamp_${STANDARDCAT}/${HEADBASE}.head  | ${P_GAWK} '{print $3}'`
      value ${CRPIX2}
      writekey ${TEMPDIR}/dummy.fits_$$ CRPIX2 "${VALUE} / WCS Coordinate reference pixel" REPLACE


      ${P_SKY2XY} -x 0.0 0.0 ${TEMPDIR}/dummy.fits_$$ @${TEMPDIR}/${BASE}_ldac.asc | ${P_GAWK} '{print $5, $6}' \
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
      ${P_LDACDELTAB} -i ${TEMPDIR}/tmp.cat1_$$\
                      -o ${TEMPDIR}/tmp.cat11_$$ -t FIELDS
      ${P_LDACRENKEY} -i ${TEMPDIR}/tmp.cat11_$$ \
                      -o ${TEMPDIR}/tmp.cat12_$$ \
                      -t STDTAB -k A_WORLD A_WCS B_WORLD B_WCS \
                                   THETA_J2000 THETAWCS
      ${P_LDACADDKEY} -i ${TEMPDIR}/tmp.cat12_$$ -t STDTAB \
                      -o ${TEMPDIR}/tmp.cat13_$$ \
		      -k CHIP ${CHIP} SHORT ""
      BADCCD=`${P_DFITS} /${MD}/${SD}/${BASE}.fits | fitsort BADCCD | grep ${BASE} | ${P_GAWK} '{print $2}'`
      ${P_LDACADDKEY} -i ${TEMPDIR}/tmp.cat13_$$ -t STDTAB \
		      -o ${MD}/${SD}/cat/${BASE}_photprep.cat \
		      -k BADCCD ${BADCCD} SHORT "Is_CCD_Bad_(1=Yes)"
    done
  else # if [ "${CATS}" != "" ]
    theli_warn "No catalogues for Chip ${CHIP} available!" "${!#}"
  fi
done

# clean up and bye
rm  ${TEMPDIR}/*_$$

theli_end "${!#}"
exit 0;
