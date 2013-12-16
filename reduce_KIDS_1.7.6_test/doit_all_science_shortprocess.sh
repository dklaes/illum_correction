#!/bin/bash -xv

# wrapper to perform all steps to process OMEGACAM science exposures
# when superflat images, reg files etc. already are present.

# you need to setup correctly at least MD, WWWDIR, SD and SURVEY below!

RUN=$1
FILT=$2

#MD=/export/euclid2_ssd/terben/pointings
MD=/vol/euclid1/euclid1_raid1/dklaes/data/
#WWWDIR=/export/euclid2_1/terben/pointings
WWWDIR=/vol/euclid1/euclid1_raid1/dklaes/data/WWW/
SD=/vol/braid1/vol3/thomas/
SURVEY=KIDS_V0.5

# set logging variables if they are not yet set:
: ${THELI_LOGGING:="Y"}
: ${THELI_LOGDIR="`pwd`/logs_${RUN}_${FILT}"}
: ${THELI_LOGLEVEL="2"}
export THELI_LOGGING
export THELI_LOGDIR
export THELI_LOGLEVEL

./doall_run_OMEGACAM_single.sh -b ${MD} \
  -s ${SURVEY} -m "UNZIPSCIENCE PREPARESCIENCE CTCOEFFSCIENCE" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/ \
  -ctd ${SD}/${SURVEY}/CT

./doall_run_OMEGACAM_single.sh -b ${MD} \
  -s ${SURVEY} -m "BIASCOPY DARKCOPY FLATCOPY SCIENCECOPY PHOTINFOCOPY" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -scb ${SD}/${SURVEY}/${FILT} -stb ${SD}/${SURVEY}/${FILT} -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "SCIENCE GLOBALWEIGHTS" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "ILLUMAPPLY" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "TESTREG WEIGHTSSCIENCE CROSSTALKWEIGHTS" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "CORRECTWEIGHTS PHOTOMSCIENCE SINGLEASTROMSCIENCE" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -sd ${SD}/${SURVEY}/

./doall_run_OMEGACAM_single.sh -b ${MD}/ \
  -s ${SURVEY} -m "SKYSUB RUNDISTRIBUTE" \
  -f ${FILT} -r ${RUN} -t SKYFLAT -www ${WWWDIR}/ \
  -bb ${SD}/${SURVEY}/BIAS -db ${SD}/${SURVEY}/DARK \
  -fb ${SD}/${SURVEY}/${FILT} \
  -sd ${SD}/${SURVEY}/ -sb ${MD}

