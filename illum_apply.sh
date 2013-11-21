#!/bin/bash

# ----------------------------------------------------------------
# File Name:           illum_apply.sh
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    31.05.2013
# Version:		V1.0
# Description:         Apply illumination correction to data
# ----------------------------------------------------------------

# Changes from V1.0 to V1.1
# - more comments

# $1  main dir
# $2  dir to apply (e.g. SCIENCE)
# $3  standard dir
# $4  filter as used in file names (e.g. r_SDSS)
# $5  extension of images
# $6  operation mode ("RUNCALIB" for illum correction for the entire run or "NIGHTCALIB" for illum correction for every night)
# $7  number of processors to be used
# $8  illum dir, where the illumination correction files are stored

MAIND=$1
APPLYD=$2
STANDARDD=$3
FILTERNAME=$4
EXTENSION=$5
MODE=$6
NPROC=$7
ILLUMDIR=$8

# including some important files
. ${INSTRUMENT:?}.ini
. ./bash_functions.include
. ./progs.ini

theli_start "$*"

if [ $# -ne 8 ]; then
  theli_error "Wrong number of command line arguments!"
  exit 1;
fi

# Check if the old not illumination corrected science image folder exist. If so, delete it and if not, create it.
if [ -d "${MAIND}/${APPLYD}/${EXTENSION}_IMAGES" ]; then
  theli_error "Illumination correction already applied! Please delete these files and move the ${EXTENSION} files back!"
  exit 1;
fi

# Check if the ILLUMDIR directory exists.
if [ ! -d "${ILLUMDIR}/${calib}/" ]; then
  theli_error "Illumination correction directory does not exist! Exiting!"
  exit 1;
fi

# Getting the nights
REDDIR=`pwd`
cd ${ILLUMDIR}/calib/
if [ "$6" == "RUNCALIB" ]; then
  NIGHTS=0
elif [ "$6" == "NIGHTCALIB" ]; then
  NIGHTS=`ls -1 | grep residuals | grep -v residuals_0 | awk -F"_" '{print $2}'`
else
  theli_error "RUNMODE not set correctly!"
  exit 1;
fi

# Check if there is at least one folder with illumination correction files...
for NIGHT in ${NIGHTS}
do
  if [ `ls ${ILLUMDIR}/calib/residuals_${NIGHT}/chip_*.fits | wc -l` -ne "${NCHIPS}" ]; then
    theli_error "No files for illumination correction avaiable!"
    exit 1;
  fi
done


cd ${REDDIR}

mkdir ${MAIND}/${APPLYD}/${EXTENSION}_IMAGES

# print out information and split the files to the available processors:
NFILES=`find ${MAIND}/${APPLYD}/ -maxdepth 1 -name \*${EXTENSION}.fits | wc -l`

echo "${NFILES} images to process with ${NPROC} processors!"

i=0
while [ $i -lt ${NPROC} ]
do
  cd ${MAIND}/${APPLYD}/
  # the 'sort' in the following pipeline ensures that the order of files is
  # always the same. This is not ensured by the 'find' command!
  FILES[$i]=`find -maxdepth 1 -name \*${EXTENSION}.fits |\
             sort | awk '(NR % '${NPROC}' == '$i')'`
  NFILES[$i]=`echo ${FILES[$i]} | wc -w`

  j=$(( $i + 1 ))

  echo -e "Starting Job ${j}. It has ${NFILES[$i]} files to process!\n"

  # launch correction job! For each processor we start a background job
  # by grouping several shell commands to a single one!
  {
    k=1
    for file in ${FILES[$i]}
    do
      BASE=`basename ${file} .fits`
      CHIP=`basename ${file} ${EXTENSION}.fits | ${P_GAWK} -F "_" '{print $NF}'`

      if [ "${MODE}" == "RUNCALIB" ]
      then
	NIGHT=0
      elif [ "${MODE}" == "NIGHTCALIB" ]
      then
	NIGHT=`${P_DFITS} ${file} | ${P_FITSORT} -d GABODSID | awk '{print $2}'`
      fi

      echo "Job ${j} (${k}/${NFILES[$i]}): processing ${file} in ${MODE} mode..."

      ${P_IC} '%1 %2 /' ${file} ${ILLUMDIR}/calib/residuals_${NIGHT}/chip_${CHIP}.fits > ${BASE}I.fits
      mv ${file} ${EXTENSION}_IMAGES/

      mv ${MAIND}/WEIGHTS/${BASE}.flag.fits ${MAIND}/WEIGHTS/${BASE}I.flag.fits
      mv ${MAIND}/WEIGHTS/${BASE}.weight.fits ${MAIND}/WEIGHTS/${BASE}I.weight.fits

      k=$(( $k + 1 ))
    done
  } &
  i=$(( $i + 1 ))
done

# only finish the script if all lauched background jobs
# are finished!
wait

theli_end
exit 0;
