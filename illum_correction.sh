#!/bin/bash

# ----------------------------------------------------------------
# File Name:           illum_correction.sh
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    02.05.2013
# Version:		V1.1
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
# + CUTRESABS with respect to MEAN

# Changes from V1.1 to V1.2
# - included _$$ to temporary files

# Changes from V1.0 to V1.1
# - corrected estimation of magnitude
# - deleted time measurements
# - deleted "illum_correction_plot_fitted.py" command because this file was combined \
#   with "illum_correction_contourplot_fitfunction.py"

MAIND=$1
STANDARDD=$2
FILTERNAME=$3
SOLUTION=$4
EXTENSION=$5
FILTER=$6
MODE=$7

MINOBJECTS=0
CUTS="RES SIGMA"
LOWERCUTPERCENT=0.1	#percent
UPPERCUTPERCENT=0.1	#percent
LOWERCUTRESABS=-0.2	#mag
UPPERCUTRESABS=0.2	#mag
LOWERCUTMAG=10		#mag
UPPERCUTMAG=25		#mag
SIGMAWIDTH=3


# including some important files
. ${INSTRUMENT:?}.ini
. ./bash_functions.include
. ./progs.ini

theli_start "$*"

# Checking for correct number of command line arguments.
if [ $# -ne 7 ]; then
  theli_error "Wrong number of command line arguments! You gave me $# but I need 7."
  exit 1;
fi

# Checking which runmode shall be used. See also information for $7 above.
REDDIR=`pwd`
cd /${MAIND}/${STANDARDD}/calib/
if [ "$7" == "RUNCALIB" ]; then
  NIGHTS=0
elif [ "$7" == "NIGHTCALIB" ]; then
  NIGHTS=`${P_LDACTOASC} -i /${MAIND}/${STANDARDD}/cat/allchips_tmp.cat -t PSSC \
      -b -k GABODSID | ${P_SORT} | uniq | awk '{printf("%s ", $1)}'`
else
  theli_error "RUNMODE not set correctly!"
  exit 1;
fi
cd ${REDDIR}

# Check, for each avaiable night, if the residuals folder exist. If so, delete it and if not, create it.
for NIGHT in ${NIGHTS}
do
  if [ ! -d "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}" ]; then
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
  else
    theli_warn "Old illumination correction detected. I will remove and redo it!"
    rm -r ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
    mkdir ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}
  fi
done

# Getting some information about the used camera. Checking first if ${CHIPGEOMETRY} is not empty.
if [ "${CHIPGEOMETRY}" == "" ]; then
  theli_error "Chipgeometry not set! Check in ${INSTRUMENT}.ini if avaiable!"
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
# and checking if they are avaiable.
CATS=`find /${MAIND}/${STANDARDD}/cat/ -name \*all_photprep_merg.cat`
if [ "${CATS}" == "" ]; then
  theli_error "No standard catalogue matched catalogues avaiable!"
  exit 1;
else
  ${P_LDACPASTE} -i ${CATS} -t PSSC\
                 -o ${TEMPDIR}/tmp_exp.cat_$$
fi

# Doing some filtering process which is the same for all chips and nights already here
# (for performance reasons).
# Filtering only those detected source with a mag in a certain filter less 99mag
${P_LDACFILTER} -i ${TEMPDIR}/tmp_exp.cat_$$ \
		  -o ${TEMPDIR}/tmp_exp.cat2_$$ \
		  -t PSSC -c "(${FILTER}mag<99);"

# Checking if we have still objects at all.
if [ ! -e "${TEMPDIR}/tmp_exp.cat2_$$" ]; then
  theli_error "No objects left after ${FILTER}mag<99 filtering!"
  exit 1;
fi


# Checking if photometric calibration file is avaiable.
if [ ! -e "${MAIND}/${STANDARDD}/calib/night_${NIGHT}_${FILTERNAME}_result.asc" ]; then
  theli_error "No photometric calibration file avaiable! Check 'night_${NIGHT}_${FILTERNAME}_result.asc'!"
  exit 1;
fi


# Now extract all needed information from the chip-based catalogues. 
# Please make sure that you have already modified stdphotom_prepare_make_ssc.conf 
# containing Xpos and Ypos (normally from input catalogue 0)!

for NIGHT in ${NIGHTS}
do
  # Getting the zeropoint with the given solution number line
  ZP=`${P_GAWK} 'NR == '${SOLUTION}' {print $1}' ${MAIND}/${STANDARDD}/calib/night_${NIGHT}_${FILTERNAME}_result.asc`
  EXT=`${P_GAWK} 'NR == '${SOLUTION}' {print $2}' ${MAIND}/${STANDARDD}/calib/night_${NIGHT}_${FILTERNAME}_result.asc`

  # First filtering for objects of this specific night. Then calculate
  # the residual with the ZP from the fit done by the THELI pipeline.
  # It's residual = detected magnitude + zeropoint - reference and
  # coordinate transformation for coordinate between -1.0 and 1.0
  if [ "${MODE}" == "RUNCALIB" ]; then
    cp ${TEMPDIR}/tmp_exp.cat2_$$ ${TEMPDIR}/tmp_exp.cat3_$$
  else
    ${P_LDACFILTER} -i ${TEMPDIR}/tmp_exp.cat2_$$ \
		    -o ${TEMPDIR}/tmp_exp.cat3_$$ \
		    -t PSSC \
		    -c "(GABODSID=${NIGHT});"
  fi

  ${P_LDACCALC} -i ${TEMPDIR}/tmp_exp.cat3_$$ \
                -o ${TEMPDIR}/tmp_exp.cat4_$$ -t PSSC \
                -c "(Mag+${ZP}+${EXT}*AIRMASS);" -n MagZP "" -k FLOAT
  ${P_LDACCALC} -i ${TEMPDIR}/tmp_exp.cat4_$$ \
                -o ${TEMPDIR}/tmp_exp.cat5_$$ -t PSSC \
                -c "(MagZP-${FILTER}mag);" -n Residual "" -k FLOAT \
                -c "((2.0*Xpos_global)/${PIXXMAX});" -n Xpos_mod "" -k FLOAT \
                -c "((2.0*Ypos_global)/${PIXYMAX});" -n Ypos_mod "" -k FLOAT


  # Now filtering according to given methods and values:
  i=0
  cp ${TEMPDIR}/tmp_exp.cat5_$$ ${TEMPDIR}/tmp_filter.cat${i}_$$
  for METHOD in ${CUTS}
  do
    i=$(( $i + 1 ))
    if [ "${METHOD}" == "PERCENT" ]; then
      # Throw away the upper and lower e.g. 10 percent (controlled via ${UPPERCUTPERCENT}
      # and ${LOWERCUTPERCENT}).
      NUMBERUPPER=`${P_LDACTOASC} -b -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC -k Residual | wc -l | ${P_GAWK} '{printf "%.0f", $1*'${UPPERCUTPERCENT}'}'`
      NUMBERLOWER=`${P_LDACTOASC} -b -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC -k Residual | wc -l | ${P_GAWK} '{printf "%.0f", $1*'${LOWERCUTPERCENT}'}'`
      LOWERVALUE=`${P_LDACTOASC} -b -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC -k Residual | sort -g | ${P_GAWK} 'NR=='${NUMBERLOWER}' {print $0}'`
      HIGHERVALUE=`${P_LDACTOASC} -b -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC -k Residual | sort -rg | ${P_GAWK} 'NR=='${NUMBERUPPER}' {print $0}'`

      ${P_LDACFILTER} -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
		      -o ${TEMPDIR}/tmp_filter.cat${i}_$$ \
		      -c "((Residual<${HIGHERVALUE})AND(Residual>${LOWERVALUE}));"


    elif [ "${METHOD}" == "RES" ]; then
      ${P_LDACFILTER} -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                      -o ${TEMPDIR}/tmp_filter.cat${i}_$$ \
                      -c "((Residual>${LOWERCUTRESABS})AND(Residual<${UPPERCUTRESABS}));"


    elif [ "${METHOD}" == "MAG" ]; then
      ${P_LDACFILTER} -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                      -o ${TEMPDIR}/tmp_filter.cat${i}_$$ \
                      -c "((MagZP>${LOWERCUTMAG})AND(MagZP<${UPPERCUTMAG}));"


    elif [ "${METHOD}" == "SIGMA" ]; then
      ${P_LDACTOASC} -b -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ -t PSSC \
                     -k Residual >> ${TEMPDIR}/res_${NIGHT}.csv_$$

      SIGMA=`${P_GAWK} '{if ($1!="#") {print $1}}' ${TEMPDIR}/res_${NIGHT}.csv_$$ \
	      | ${P_GAWK} -f meanvar.awk | grep sigma | ${P_GAWK} '{print $3}'`
      MEAN=`${P_GAWK} '{if ($1!="#") {print $1}}' ${TEMPDIR}/res_${NIGHT}.csv_$$ \
	      | ${P_GAWK} -f meanvar.awk | grep mean | ${P_GAWK} '{print $3}'`

      # Filtering only those residuals which lies in certain given limits.
      ${P_LDACFILTER} -i ${TEMPDIR}/tmp_filter.cat$(( $i - 1 ))_$$ \
                      -o ${TEMPDIR}/tmp_filter.cat${i}_$$ \
                      -c "((Residual<${MEAN}+${SIGMAWIDTH}*${SIGMA})AND(Residual>${MEAN}-${SIGMAWIDTH}*${SIGMA}));" \
                      -t PSSC
    fi
  done

  cp ${TEMPDIR}/tmp_filter.cat${i}_$$ ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat

  # Splitting up one catalogue with all chips into ${NUMCHIPS} files.
  # Check, if for all chips enough objects are avaiable. If not, abort.
  # Extracting all needed information into a CSV file (night based)
  i=1
  while [ ${i} -le ${NCHIPS} ]
  do
    ${P_LDACFILTER} -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all_filtered.cat -t PSSC \
                    -o ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered.cat \
                    -c "(CHIP=${i});"
    
    if [ -e "${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered.cat" ]; then
      NUMBER=`${P_LDACTOASC} -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered.cat -t PSSC -k MagZP | wc -l`
      if [ ${NUMBER} -le ${MINOBJECTS} ]; then
	theli_error "Not enough objects avaiable for fitting. Chip ${CHIP} caused the first problem!"
	exit 1;
      fi
    else
      theli_error "No information for at least one chip avaiable. Chip ${CHIP} caused the first problem!"
      exit 1;
    fi

    ${P_LDACTOASC} -b -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered.cat -t PSSC \
	    -k Xpos Ypos Mag MagErr ${FILTER}mag IMAGEID Residual Xpos_mod \
	    Ypos_mod AIRMASS Xpos_global Ypos_global MagZP >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv
    i=$(( $i + 1 ))
  done

  
  # Fitting the data
  ./illum_correction_fit.py ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/

  rm ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_*.csv
  
  #Applying the fit-parameter to our catalog data...
  i=1
  # Getting the prefaktors of our correction model.
  A=`grep A ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
  B=`grep B ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
  C=`grep C ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
  D=`grep D ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
  E=`grep E ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
  while [ ${i} -le ${NCHIPS} ]
  do
    FCHIP=`grep -m1 F${i} ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat | ${P_GAWK} '{print $3}'`
    
    # Calculating the fitted magnitude for each object and the fitted residuals.
    ${P_LDACCALC} -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered.cat -o ${TEMPDIR}/chip_${i}_filtered_fitted.cat1_$$ -t PSSC \
      -c "(MagZP-${A}*(Xpos_mod*Xpos_mod)-${B}*(Ypos_mod*Ypos_mod)-${C}*(Xpos_mod*Ypos_mod)-${D}*Xpos_mod-${E}*Ypos_mod-${FCHIP});" \
      -n Mag_fitted "" -k FLOAT
    ${P_LDACCALC} -i ${TEMPDIR}/chip_${i}_filtered_fitted.cat1_$$ -o ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered_fitted.cat -t PSSC \
      -c "(Mag_fitted-${FILTER}mag);" -n Residual_fitted "" -k FLOAT
    
    # Extracting all needed information into a CSV file (night based)
    ${P_LDACTOASC} -b -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_filtered_fitted.cat -t PSSC -k Xpos Ypos Mag MagErr ${FILTER}mag IMAGEID \
      Residual Xpos_mod Ypos_mod Mag_fitted Residual_fitted AIRMASS Xpos_global Ypos_global MagZP >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv
    i=$(( $i + 1 ))
  done

  # Now calculating some statistics for before and after the fitting process for each chip and for all chips...
  i=1
  while [ ${i} -le ${NCHIPS} ]
  do
    echo "" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    echo "Statistics of residuals before fitting:" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    ${P_GAWK} '{if ($1!="#") {print $7}}' ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv | ${P_GAWK} -f meanvar.awk >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    echo "" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    echo "Statistics of residuals after fitting:" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    ${P_GAWK} '{if ($1!="#") {print $11}}' ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv | ${P_GAWK} -f meanvar.awk >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.dat
    cat ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.csv
    i=$(( $i + 1 ))
  done
  echo "" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat
  echo "Statistics of residuals before fitting:" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat
  ${P_GAWK} '{if ($1!="#") {print $7}}' ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.csv | ${P_GAWK} -f meanvar.awk >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat
  echo "" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat
  echo "Statistics of residuals after fitting:" >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat
  ${P_GAWK} '{if ($1!="#") {print $11}}' ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.csv | ${P_GAWK} -f meanvar.awk >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.dat &
  
  # Creating a contour plot from the correction function and create a correction FITS file...
  SIGMA=`${P_GAWK} '{if ($1!="#") {print $7}}' ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_all.csv \
    | ${P_GAWK} -f meanvar.awk | grep sigma | ${P_GAWK} '{print $3}'`
  ./illum_correction_contourplot_fitfunction.py ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/ -${SIGMA} ${SIGMA}
done


# Cleaning up...
rm ${TEMPDIR}/*_$$ &

theli_end
exit 0;
