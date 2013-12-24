SIGMA=0.0503426
MAIND=/home/dklaes/git/illum_correction/testing/run_11_08_f/
STANDARDD=STANDARD_r_SDSS
NIGHT=0

. ./OMEGACAM.ini

./illum_correction_contourplot_fitfunction.py ${MAIND}/${STANDARDD}/calib/residuals_${NIGHT}/ -${SIGMA} ${SIGMA}
