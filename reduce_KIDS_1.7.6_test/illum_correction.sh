#!/bin/bash

# ----------------------------------------------------------------
# File Name:           illum_correction.sh
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    26.11.2013
# Version:	       V1.2
# Description:         Estimate position dependent magnitude offset
# ----------------------------------------------------------------


# $1			main dir
# $2			standard dir
# $3			filter as used in file names (e.g. r_SDSS)
# $4			Solution (row) number from calib/....asc
# $5			extension of images
# $6			filter as used in files (e.g. R instead of r_SDSS)
# $7  			operation mode ("RUNCALIB" for illum correction for the entire
#     			run or "NIGHTCALIB" for illum correction for every night)
# $8			name of color term (e.g. gmr)
#
# MINOBJECTS		Minimal number of objects that are required for fitting
# CUTS			Array containing cut types and order
#			At the moment: PERCENT, RES, MAG and SIGMA
# LOWERCUTPERCENT	The lower part of the data is thrown away right in the beginning (given in percent)
# UPPERCUTPERCENT	The upper part of the data is thrown away right in the beginning (given in percent)
# LOWERCUTRESABS	Residuals with smaller residuals than this value will be cutted
# UPPERCUTRESABS	Residuals with larger residuals than this value will be cutted
# LOWERCUTMAG		Minimal magnitude being considered
# UPPERCUTMAG		Maximal magnitude being considered
# SIGMAWIDTH		How many sigmas shall be taken?
# LOWERCUTRESMEAN	Residuals with smaller residuals than this value with respect to mean will be cutted
# UPPERCUTRESMEAN	Residuals with larger residuals than this value with respect to mean will be cutted

# Changes from V1.0 to V1.1
# - corrected estimation of magnitude
# - deleted time measurements
# - deleted "illum_correction_plot_fitted.py" command because this file was combined \
#   with "illum_correction_contourplot_fitfunction.py"

# Changes from V1.1 to V1.2
# - included _$$ to temporary files

# Changes from V1.2 to V1.3
# - included residual cut with respect to the mean
# - Single plots are available
# - Added color term name and include color term
# - Included error calculation for residual with Gaussian error propagation


MAIND=$1
STANDARDD=$2
FILTERNAME=$3
SOLUTION=$4
EXTENSION=$5
FILTER=$6
MODE=$7
COLOR=$8

MINOBJECTS=0
CUTS="RESMEAN"
LOWERCUTPERCENT=10	#percent
UPPERCUTPERCENT=10	#percent
LOWERCUTRESABS=-0.2	#mag
UPPERCUTRESABS=0.2	#mag
LOWERCUTMAG=10		#mag
UPPERCUTMAG=25		#mag
SIGMAWIDTH=1
LOWERCUTRESMEAN=-0.2	#mag
UPPERCUTRESMEAN=0.2	#mag


# including some important files
. ${INSTRUMENT:?}.ini
. ./bash_functions.include
. ./progs.ini

theli_start "$*"

# Checking for correct number of command line arguments.
if [ $# -ne 8 ]; then
  theli_error "Wrong number of command line arguments! You gave me $# but I need 8."
  exit 1;
fi

# Checking which runmode shall be used. See also information for $7 above.
if [ "$7" == "RUNCALIB" ]; then
  NIGHTS=0
elif [ "$7" == "NIGHTCALIB" ]; then
  NIGHTS=`${P_PYTHON} illum_ldactools.py -i /${MAIND}/${STANDARDD}/cat/allexp_tmp.cat -t PSSC \
	-k GABODSID -a UNIQUE_ELEMENTS`
else
  theli_error "RUNMODE not set correctly!"
  exit 1;
fi


# Check, for each available night, if the residuals folder exist. If so, delete it and if not, create it.
for NIGHT in ${NIGHTS}
do
  if [ ! -d "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}" ]; then
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/single_plots
  else
    theli_warn "Old illumination correction detected. I will remove and redo it!"
    rm -r ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/single_plots
  fi
done

# Getting some information about the used camera. Checking first if ${CHIPGEOMETRY} is not empty.
if [ "${CHIPGEOMETRY}" == "" ]; then
  theli_error "Chipgeometry not set! Check in ${INSTRUMENT}.ini if available!"
  exit 1;
else
  ROWMAX=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $1}'`
  COLUMNMAX=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $2}'`
  PIXX=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $3}'`
  PIXY=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $4}'`
  PIXXMAX=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $1*$3}'`
  PIXYMAX=`echo ${CHIPGEOMETRY} | ${P_GAWK} '{print $2*$4}'`
fi

# Getting all cats into one catalogue (modification needed to update from reduce_KIDS_1.7.2 to 1.7.6)
# and checking if they are available.
CATS=`find /${MAIND}/${STANDARDD}/cat/ -name \*all_photprep_merg.cat`
if [ "${CATS}" == "" ]; then
theli_error "No standard catalogue matched catalogues available!"
  exit 1;
else
  ${P_PYTHON} illum_ldactools.py -i "${CATS}" -t PSSC -a PASTE_TABLES \
                 -o ${TEMPDIR}/tmp_exp.cat_$$
fi

# Doing some filtering process which is the same for all chips and nights already here
# (for performance reasons).
# Filtering only those detected source with a mag in a certain filter less 99mag
${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp.cat_$$ \
		  -o ${TEMPDIR}/tmp_exp_filter.cat_$$ \
		  -t PSSC -a FILTER_USUABLE -e "${FILTER} ${COLOR}"

# Checking if we have still objects at all.
if [ ! -e "${TEMPDIR}/tmp_exp_filter.cat_$$" ]; then
  theli_error "No temporary file after ${FILTER}mag<99 filtering available!"
  exit 1;
else
  NUMOBJECTS=`${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp_filter.cat_$$ -t PSSC -a NUMBER_OF_ELEMENTS`
  if [ ${NUMOBJECTS} -eq 0 ]; then
	  theli_error "No objects left after ${FILTER}mag<99 filtering!"
	  exit 1;
  fi
fi

${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp_filter.cat_$$ \
		  -o ${TEMPDIR}/tmp_exp.cat2_$$ \
		  -a FILTER_ELEMENTS \
		  -t PSSC -k BADCCD \
		  -c "=" -v 0

# Checking if we have still objects at all.
if [ ! -e "${TEMPDIR}/tmp_exp.cat2_$$" ]; then
  theli_error "No temporary file after BADCCD flag filtering available!"
  exit 1;
else
  NUMOBJECTS=`${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp.cat2_$$ -t PSSC -a NUMBER_OF_ELEMENTS`
  if [ ${NUMOBJECTS} -eq 0 ]; then
	  theli_error "No objects left after BADCCD flag filtering!"
	  exit 1;
  fi
fi


# Checking if photometric calibration file is available.
for NIGHT in ${NIGHTS}
do
	if [ ! -e "${MAIND}/${STANDARDD}/calib/night_${NIGHT}_${FILTERNAME}_result.asc" ]; then
	  theli_error "No photometric calibration file available! Check 'night_${NIGHT}_${FILTERNAME}_result.asc'!"
	  exit 1;
	fi
done


# Now extract all needed information from the chip-based catalogues.
# Please make sure that you have already modified stdphotom_prepare_make_ssc.conf
# containing Xpos and Ypos (normally from input catalogue 0)!

for NIGHT in ${NIGHTS}
do
  # First filtering for objects of this specific night. Then calculate
  # the residual with the ZP from the fit done by the THELI pipeline.
  # It's residual = detected magnitude + zeropoint - reference and
  # coordinate transformation for coordinate between -1.0 and 1.0
  if [ "${MODE}" == "RUNCALIB" ]; then
    cp ${TEMPDIR}/tmp_exp.cat2_$$ ${TEMPDIR}/tmp_exp.cat3_$$
  else
    ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp.cat2_$$ \
		    -o ${TEMPDIR}/tmp_exp.cat3_$$ \
		    -a FILTER_ELEMENTS \
		    -t PSSC -k GABODSID \
		    -c "=" -v ${NIGHT}
  fi

  ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_exp.cat3_$$ \
			-o ${TEMPDIR}/tmp_exp.cat5_$$ -t PSSC \
			-a CALCS_BEFORE_FITTING \
			-e "${MAIND}/${STANDARDD}/calib/night_${NIGHT}_${FILTERNAME}_result.asc 2 ${COLOR} ${FILTER}"

  # Now filtering according to given methods and values:
  i=0
  cp ${TEMPDIR}/tmp_exp.cat5_$$ ${TEMPDIR}/tmp_filter.cat${i}_$$
  for METHOD in ${CUTS}
  do
    i=$(( $i + 1 ))
    if [ "${METHOD}" == "PERCENT" ]; then
      # Throw away the upper and lower e.g. 10 percent (controlled via ${UPPERCUTPERCENT}
      # and ${LOWERCUTPERCENT}).
      ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
		      -o ${TEMPDIR}/tmp_filter.cat${i}_$$ -a FILTER_PERCENT \
		      -e "${LOWERCUTPERCENT} ${UPPERCUTPERCENT}" -k Residual


    elif [ "${METHOD}" == "RES" ]; then
	${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                       -o ${TEMPDIR}/tmp_filter.cat${i}_$$ -a FILTER_RESIDUAL \
		        -e "${LOWERCUTRESABS} ${UPPERCUTRESABS}" -k Residual


    elif  [ "${METHOD}" == "RESMEAN" ]; then
      ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                       -o ${TEMPDIR}/tmp_filter.cat${i}_$$ -a FILTER_RESIDUALMEAN \
			-e "${LOWERCUTRESMEAN} ${UPPERCUTRESMEAN}" -k Residual


    elif [ "${METHOD}" == "MAG" ]; then
      ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                       -o ${TEMPDIR}/tmp_filter.cat${i}_$$ -a FILTER_MAGNITUDE \
			-e "${LOWERCUTMAG} ${UPPERCUTMAG}" -k MagZP


    elif [ "${METHOD}" == "SIGMA" ]; then
      ${P_PYTHON} illum_ldactools.py -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                       -o ${TEMPDIR}/tmp_filter.cat${i}_$$ -a FILTER_SIGMA \
			-e ${SIGMAWIDTH} -k Residual
    fi
  done

  cp ${TEMPDIR}/tmp_filter.cat${i}_$$ ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat

  # Splitting up one catalogue with all chips into ${NUMCHIPS} files.
  # Check, if for all chips enough objects are available. If not, warn.
  i=1
  if [ -e "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat" ]; then
	  while [ ${i} -le ${NCHIPS} ]
	  do
	      NUMBER=`${P_PYTHON} illum_ldactools.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat -t PSSC -a CHECK_ENOUGH_OBJECTS -v ${i}`

	      if [ ${NUMBER} -le ${MINOBJECTS} ]; then
		theli_warn "Not enough objects available for fitting. Chip ${i} caused the problem!"
	      fi
              i=$(( $i + 1 ))
	  done
  else
	theli_warn "No information for all chips available!"
  fi

  # Fitting the data
  ${P_PYTHON} illum_correction_fit.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat \
				-t PSSC -p ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/

  #Applying the fit-parameter to our catalog data...

  ${P_PYTHON} illum_ldactools.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat \
			-o ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered_fitted.cat -t PSSC \
			-a CALCS_AFTER_FITTING \
			-e "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/coeffs.txt ${FILTER}"

  ${P_PYTHON} illum_ldactools.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered_fitted.cat -t PSSC \
			-a STATISTICS -e "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/coeffs.txt 10 10" \
			-o ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/stats.txt

  ${P_PYTHON} illum_ldactools.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered_fitted.cat -t PSSC \
			-o ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/fitting.txt -a MAG_DEPENDENCY

  # Creating a contour plot from the correction function and create a correction FITS file...
  ${P_PYTHON} illum_correction_contourplot_fitfunction.py -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered_fitted.cat -t PSSC \
			-p ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/ -e "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/coeffs.txt"
done


# Cleaning up...
rm ${TEMPDIR}/*_$$

theli_end
exit 0;
