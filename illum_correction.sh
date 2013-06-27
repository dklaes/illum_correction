#!/bin/bash

# ----------------------------------------------------------------
# File Name:           illum_correction.sh
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    02.05.2013
# Version:		V1.1
# Description:         Estimate position dependent magnitude offset
# ----------------------------------------------------------------


# $1  main dir
# $2  standard dir
# $3  filter as used in file names (e.g. r_SDSS)
# $4  Solution (row) number from calib/....asc
# $5  extension of images
# $6  filter as used in files (e.g. R instead of r_SDSS)
# $7  operation mode ("RUNCALIB" for illum correction for the entire 
#     run or "NIGHTCALIB" for illum correction for every night)


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

# Splitting up one catalogue with all chips into ${NUMCHIPS} files.
i=1
  while [ ${i} -le ${NCHIPS} ]
  do
	${P_LDACFILTER} -i ${TEMPDIR}/tmp_exp.cat2_$$ -t PSSC \
			-o ${MAIND}/${STANDARDD}/cat/chip_${i}_merg.cat \
			-c "(CHIP=${i});"
  i=$(( $i + 1 ))
done

# Check, if for all chips are information avaiable. If not, abort.
for CHIP in ${NUMCHIPS}
do
  if [ -e "${MAIND}/${STANDARDD}/cat/chip_${i}_merg.cat" ]; then
    NUMBER=`${P_LDACDESC} -i ${MAIND}/${STANDARDD}/cat/chip_${i}_merg.cat | grep elements | ${P_GAWK} -F . '{print $1}'`
    if [ ${NUMBER} -le 0 ]; then
      theli_error "Not enough objects avaiable for fitting. Chip ${CHIP} caused the first problem!"
      exit 1;
    fi
  else
    theli_error "No information for all chips avaiable. Chip ${CHIP} caused the first problem!"
    exit 1;
  fi
done

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
  i=1
  while [ ${i} -le ${NCHIPS} ]
  do
    # Calculating the residual with the ZP from the fit done by the THELI pipeline.
    # It's residual = detected magnitude + zeropoint - reference and
    # coordinate transformation for coordinate between -1.0 and 1.0
    ${P_LDACCALC} -i ${MAIND}/${STANDARDD}/cat/chip_${i}_merg.cat \
	    -o ${TEMPDIR}/chip_${i}_merg_corr1.cat_$$ -t PSSC \
	    -c "(Mag+${ZP}+${EXT}*AIRMASS);" -n MagZP "" -k FLOAT
    ${P_LDACCALC} -i ${TEMPDIR}/chip_${i}_merg_corr1.cat_$$ \
	    -o ${TEMPDIR}/chip_${i}_merg_corr2.cat_$$ -t PSSC \
	    -c "(MagZP-${FILTER}mag);" -n Residual "" -k FLOAT \
	    -c "((2.0*Xpos_global)/${PIXXMAX});" -n Xpos_mod "" -k FLOAT \
	    -c "((2.0*Ypos_global)/${PIXYMAX});" -n Ypos_mod "" -k FLOAT
    ${P_LDACTOASC} -b -i ${TEMPDIR}/chip_${i}_merg_corr2.cat_$$ -t PSSC \
	    -k Residual >> ${TEMPDIR}/res_${NIGHT}.csv_$$
    i=$(( $i + 1 ))
  done

  SIGMA=`${P_GAWK} '{if ($1!="#") {print $1}}' ${TEMPDIR}/res_${NIGHT}.csv_$$ \
	  | ${P_GAWK} -f meanvar.awk | grep sigma | ${P_GAWK} '{print $3}'`
  MEAN=`${P_GAWK} '{if ($1!="#") {print $1}}' ${TEMPDIR}/res_${NIGHT}.csv_$$ \
	  | ${P_GAWK} -f meanvar.awk | grep mean | ${P_GAWK} '{print $3}'`
  

  i=1
  while [ ${i} -le ${NCHIPS} ]
  do
    # Filtering only those residuals which lies in certain given limits.
    ${P_LDACFILTER} -i ${TEMPDIR}/chip_${i}_merg_corr2.cat_$$ \
	    -o ${TEMPDIR}/chip_${i}_merg_corr.cat_$$ -t PSSC \
	    -c "((Residual<${MEAN}+${SIGMA})AND(Residual>${MEAN}-${SIGMA}));"
    # Extracting all needed information into a CSV file (night based)
    ${P_LDACTOASC} -b -i ${TEMPDIR}/chip_${i}_merg_corr.cat_$$ -t PSSC \
	    -k Xpos Ypos Mag MagErr ${FILTER}mag IMAGEID Residual Xpos_mod \
	    Ypos_mod AIRMASS Xpos_global Ypos_global MagZP >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv
    mv ${TEMPDIR}/chip_${i}_merg_corr.cat_$$ ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_merg_corr.cat &
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
    ${P_LDACCALC} -i ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_merg_corr.cat -o ${TEMPDIR}/chip_${i}_merg_fitted1.cat_$$ -t PSSC \
      -c "(MagZP-${A}*(Xpos_mod*Xpos_mod)-${B}*(Ypos_mod*Ypos_mod)-${C}*(Xpos_mod*Ypos_mod)-${D}*Xpos_mod-${E}*Ypos_mod-${FCHIP});" \
      -n Mag_fitted "" -k FLOAT
    ${P_LDACCALC} -i ${TEMPDIR}/chip_${i}_merg_fitted1.cat_$$ -o ${TEMPDIR}/chip_${i}_merg_fitted.cat_$$ -t PSSC \
      -c "(Mag_fitted-${FILTER}mag);" -n Residual_fitted "" -k FLOAT
    
    # Extracting all needed information into a CSV file (night based)
    ${P_LDACTOASC} -b -i ${TEMPDIR}/chip_${i}_merg_fitted.cat_$$ -t PSSC -k Xpos Ypos Mag MagErr ${FILTER}mag IMAGEID \
      Residual Xpos_mod Ypos_mod Mag_fitted Residual_fitted AIRMASS Xpos_global Ypos_global MagZP >> ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}.csv
    mv ${TEMPDIR}/chip_${i}_merg_fitted.cat_$$ ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/chip_${i}_merg_fitted.cat &
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
  ./illum_correction_contourplot_fitfunction.py ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/ -${SIGMA} ${SIGMA}
done


# Cleaning up...
rm ${TEMPDIR}/*_$$ &

theli_end
exit 0;
