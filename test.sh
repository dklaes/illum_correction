#!/bin/bash -xv

SIGMA=0.0503426
MAIND=/home/dklaes/git/illum_correction/testing/run_11_08_f/
STANDARDD=STANDARD_r_SDSS
NIGHT=0
NPARA=1

. ./OMEGACAM.ini

rm testing/run_11_08_f/STANDARD_r_SDSS/calib/residuals_0/*.fits
rm testing/run_11_08_f/STANDARD_r_SDSS/calib/residuals_0/*.png
rm testing/run_11_08_f/STANDARD_r_SDSS/calib/residuals_0/single_plots/*.png

python illum_correction_contourplot_fitfunction.py ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/ -${SIGMA} ${SIGMA}
