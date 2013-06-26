This version V1.1 was tested with ldacpipeline-1.7.6 and reduce_KIDS_1.7.6

CHANGES to V1.0:
- create_stdphotom_prepare_para.sh:
-- Added CHIP value to catalog
-- Modified it to the new version of file

- illum_apply.sh:
-- More comments

- illum_correction.sh:
-- corrected estimation of magnitude with colour term
-- deleted time measurements
-- deleted "illum_correction_plot_fitted.py" command because this file was combined \
   with "illum_correction_contourplot_fitfunction.py"

- illum_correction_contourplot_fitfunction.py:
-- more comments
-- Checkplot magnitude / residual dependency added
-- Now running also without X-server
-- old "illum_correction_plot_fitted.py" included and restructed to multiprocessing

- illum_correction_fit.py
-- more comments



CHANGES to be done to make it work:

reduce_KIDS_1.7.6:
- progs.ini:
-- Add "P_SKY2XY=..." for sky to xy coordinate transformation
-- Add "P_XY2SKY=..." for xy to sky coordinate transformation
-- Modify "P_ACLIENT=..." to find correct program
-- Modify "P_SM=..." to find correct program
-- Add "export NPARA" so it can be used by other programs in the terminal

- create_stdphotom_prepare_para.sh: Replace
-- Adding "Xpos", "Ypos", "Xpos_global", "Ypos_global" and "CHIP" from the first catalog

- create_stdphotom_merge_exposurepara.sh
-- Modified for global coordinates

- doall_run_OMEGACAM_single.sh: Replace
-- Adding commands for illumcorrection and illumapply

- meanvar.awk: copy
-- Calculate mean, cariance, standard deviation and number of datapoints as (g)awk script

- OMEGACAM.ini:
-- Add "export CHIPGEOMETRY" so it can be used by other programs in the terminal
-- Add comment (see [1]) for better understanding
-- Add "OFFSETX=..." (for OMEGACAM use [2]), see also [1]
-- Add "OFFSETY=..." (for OMEGACAM use [3]), see also [1]
-- Add "export OFFSETX" so it can be used by other programs in the terminal
-- Add "export OFFSETY" so it can be used by other programs in the terminal


ldacpipeline-1.7.6:
- stdphotom_prepare_global_coordinates.conf: Add to ${DATACONF}
-- Addition needed for global coordinates



[1]:	# offset needed for correct plotting and fitting of residuals
	# It is the offset of the (0,0) coordinate of the corresponding chip
	# compared to the (0,0) coordinate of the camera.
[2]:	OFFSETX="-8552 -6412 -4269 -2130 10 2152 4294 6432 -8552 -6412 -4269 -2130 10 2152 4294 6432 -8552 -6412 -4269 -2130 10 2152 4294 6432 -8552 -6412 -4269 -2130 10 2152 4294 6432"
[3]:	OFFSETY="-8583 -8583 -8583 -8583 -8583 -8583 -8583 -8583 -4124 -4124 -4124 -4124 -4124 -4124 -4124 -4124 -81 -81 -81 -81 -81 -81 -81 -81 4380 4380 4380 4380 4380 4380 4380 4380"
