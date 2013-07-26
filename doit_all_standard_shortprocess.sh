#!/bin/bash -xv

# wrapper to perform all steps to process OMEGACAM standard exposures
# when superflat images, reg files etc. already are present.

# you need to setup correctly at least MD, WWWDIR and SD below!

RUN=$1
FILT=$2

MD=/vol/euclid1/euclid1_raid1/dklaes/data/
WWWDIR=/vol/euclid1/euclid1_raid1/dklaes/data/WWW/
SD=/vol/braid1/vol3/thomas/
SURVEY=KIDS_V0.5

# set logging variables if they are not yet set:
: ${THELI_LOGGING:="Y"}
: ${THELI_LOGDIR:="`pwd`/logs_STANDARD_${RUN}_${FILT}"}
: ${THELI_LOGLEVEL:="1"}
export THELI_LOGGING
export THELI_LOGDIR
export THELI_LOGLEVEL

#./doall_run_OMEGACAM_single.sh -b ${MD} \
#  -s ${SURVEY} -m "UNZIPSTANDARD PREPARESTANDARD" \
#  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
#  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
#  -fb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/


#./doall_run_OMEGACAM_single.sh -b ${MD} \
#  -s ${SURVEY} -m "BIASCOPY DARKCOPY FLATCOPY SCIENCECOPY" \
#  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
#  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
#  -fb ${SD}/${SURVEY}/${FILT} -scb ${SD}/${SURVEY}/${FILT} \
#  -stb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/

#./doall_run_OMEGACAM_single.sh -b ${MD}/ \
#  -s ${SURVEY} -m "STANDARD GLOBALWEIGHTS ABSPHOTOM" \
#  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
#  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
#  -fb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "ABSPHOTOM" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "ILLUMCORRECTION" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/ \
  -id ${MD}/${SURVEY}/${FILT}/${RUN}/STANDARD_${FILT}/calib/


#./copy_KIDS_products.sh ${MD}/${SURVEY} ${SD}/${SURVEY} ${RUN} ${FILT} STANDARD


#./doall_run_OMEGACAM_single.sh -b ${MD}/ \
#  -s ${SURVEY} -m "CLEANSTANDARD CLEANSCIENCE CLEANCALIB DELETEDIRS" \
#  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
#  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
#  -fb ${SD}/${SURVEY}/${FILT} \
#  -sd ${SD}/${SURVEY}/

