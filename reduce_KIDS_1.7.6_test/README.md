Illumination correction for reduce_KIDS_1.7.6_test
==================================================

This repository contains all files needed for an illumination correction in the framework of THELI for reduce_KIDS_1.7.6_test.

Now dependencies:
- sky2xy
- xy2sky
- scipy-0.11
- numpy-1.6.2
- astropy-0.2.3

The following files have to be saved into the reduce directory:

- meanvar.awk
- OMEGACAM.ini (now contains the offsets between the individual chips)
- progs.ini (now contains the paths to sky2xy and xy2sky)
- illum_correction_contourplot_fitfunction.py
- illum_correction_fit.py
- create_stdphotom_merge_exposurepara.sh
- create_stdphotom_prepare_para.sh
- doall_run_OMEGACAM_single.sh
- doit_all_science_shortprocess.sh
- doit_all_standard_shortprocess.sh
- illum_apply.sh
- illum_correction.sh


The following files have to be saved into the CONF directory (path to CONF directory see progs.ini):
- stdphotom_prepare_global_coordinates.conf


Please do not forget to change all corresponding paths, especially in progs.ini and all doall* or doit_all* files!
