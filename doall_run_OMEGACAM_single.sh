#!/bin/bash -xv

# the script gives the reduction steps for a OMEGACAM run reduction;
# ldacpipeline version 1.5.3 and higher. The current version deals
# with SCIENCE, SCIENCESHORT and STANDARDSTAR observaions. The
# calibration data BIAS, DARK and FLAT are processed by other scripts.

# HISTORY:
# 07.11.2011:
# project started
#
# 06.02.2012:
# - adapation of SCIENCE image limits for influencing the global weight.
#   OMEGACAM has very low counts in u-band which leads to larger errors
#   in science images.
# - fixing of left-overs of copy/paste from CFHTLS scripts.
# - splitting up of prepare step for different image types (we often do not
#   want to do the time intensive splitting of short science frames etc.)
# - In the SINGLEASTROM step I now delete an old cat subdirectory if present.
#   I saw in previous reductions that catalogues merged from SExtractor and 
#   KSB were messed up when old catalogues were present. This, up to now,
#   did not happen with a clean catalogue creation.
#
# 28.02.2012:
# - I increased the limit for DARK values entering the globalweight/flag.
#   I saw that chip 5 can have a broad, high-valued middle vertical stripe.
# - I moved creation of normalised flats from the FLAT to the GLOBALWEIGHT
#   mode. This ensures redoing the images of gloablweights are redone.
#
# 19.11.2012:
# I started to split the functionality of this script. It now is
# only responsible for processing SCIENCE, SCIENCESHORT and STANDARD
# images. BIAS, DARK and FLATS are processed with other scripts.
#
# 18.12.2012:
# I included the possibility to process SCIENCESHORT exposures similar
# to SCIENCE data (Weight images, singleastrom catalogues). These processing
# steps are requried to use the SCIENCESHORT exposures for astrometric
# calibration later.
#
# 30.12.2012:
# The script 'pipelog.py' has a new command line argument to indicate
# photometric calibration on RUN basis.
#
# 10.01.2013:
# The script create_stdphotom_prepare.sh is now called
# create_stdphotom_prepare_para.sh because it is a parallel script!
#
# 16.01.2013:
# - I corrected a bug in the PHOTOMSCIENCESHORT mode. The zeropoints
#   were not transfered because we checked for the existence of a wrong file.
# - I introduced new modes to be able to process SCIENCESHORT images
#   independently of the SCIENCE images (but of course with necessary
#   products being available in a storage directory).
# - By default, temporary directories with not needed data are deleted
#   during processing. If this should not be done the command line
#   option '-noclean' should be used.
#
# 17.01.2013:
# The 'WEIGHTS' mode (get weights for SCIENCE and SCIENCESHORT
# exposures) also matched 'GLOBALWEIGHTS' - fixed.
#
# 18.01.2013:
# I had to substitute some 'rm' by 'rm -rf' to delete directories with
# their contents. In the GLOBALWEIGHTS mode an 'old' WEIGHTS directory
# is now deleted. This allows us to use better the diskspace during
# SCIENCE and STANDARDSTAR processing (STANDARDSTAR files can be
# removed completely once zeropoints are determined.
#
# 05.02.2013:
# I made some modifications to (re)process SCIENCE images. For KIDS
# we make a first pass thorugh science images and only store absolutely
# necessary products to be able to (re)process them quickly. The
# script was adapted for this task.
#
# 11.02.2013:
# new mode CROSSTALKWEIGHTS to correct science image weights from chip 25
# for crosstalk saturation effects imposed by chip 26.
#
# 19.02.2013:
# The superflat threshold was not set correctly in standardstar/scienceshort
# image processing - fixed.
#
# 14.03.2013:
# the masking of saturated crosstalk pixels was expanded to chips 25, 26
# and 27.
#
# 03.04.2013:
# I added crosstalk correction for science images in the PREPARE step.

# the script is called with
# ./doall_run_OMEGACAM_single.sh -m "MODE1 MODE2 ...."
# the different modes are the following:
#
# - UNZIP....
#   unzip RICE compressed original data
#
# - PREPARE......
#   splits data and checks whether the modes of images are within
#   given limits
#
# - SCIENCE
#   does the preprocessing on science frames (overscan, flatfielding,
#   superflatfielding)
#
# - SCIENCESHORT
#
# - STANDARD
#
# - CLEAN.....
#   deletes all files to bring directories to their original state
#
# - CLEANSTANDARDDIRS
#   deletes all necessary, temporaray files
#   created during the standardstar frame processing
#
# - MOSAICS
#   creates small, binned mosaics from SCIENCE and SCIENCESHORT
#   exposures
#
# - SINGLEASTROM...........  
#   perform singlastrom on run basis (this is mainly for the
#   completeness of the WEB pages to be created on a run basis)
#
# - WEIGHTS.........
#   creates weight images
#
# - ABSPHOTOM
#   perform absolute photometry on standardstar observations
#
# - RUNDISTRIBUTE
#   distribute current OMEGACAM run into the OMEGACAM/KIDS set structure
#

#
# definitions and standard values:
#
export INSTRUMENT=OMEGACAM
BASEDIR=/raid4/thomas/SABINE/I_EIS
RUN=run_2004_03_14

# band to be processed (mainly needed to decide whether fringing has
# to be done or not (filters are i_SDSS, etc)
FILTER="g_SDSS"

#
# Flatfield to be used (DOMEFLAT or SKYFLAT)
FLATTYPE=SKYFLAT

#
# by default we do superflatfielding:
SUPERFLAT=S

#
# default survey:
SURVEY=KIDS

#
# read command line arguments
#
MODE=""
NOTRACKS=0
TRACKS=0
BIASBASEDIR=""
DARKBASEDIR=""
WEBDIRROOT=""
SUPERFLLISTDIR="."
CTLISTDIR="."
SCIENCEBASEDIR=""
STANDARDBASEDIR=""
SETBASEDIR=""
CLEANUP=1

GoOn=0
while [ $GoOn = 0 ]
do
   case $1 in
   -b)
       BASEDIR=${2}
       shift 2
       ;;
   -bb)
       BIASBASEDIR=${2}
       shift 2
       ;;
   -db)
       DARKBASEDIR=${2}
       shift 2
       ;;
   -ctd)
       CTLISTDIR=${2}
       shift 2
       ;;
   -f)
       FILTER=${2}
       shift 2
       ;;
   -fb)
       FLATBASEDIR=${2}
       shift 2
       ;;
   -m)
       MODE=${2}
       shift 2
       ;;
   -noclean)
       CLEANUP=0  
       shift 1
       ;;
   -notracks)
       NOTRACKS=1  
       shift 1
       ;;
   -r)
       RUN=${2}
       shift 2
       ;;
   -s)
       SURVEY=${2}
       shift 2
       ;;
   -sb)
       SETBASEDIR=${2}
       shift 2
       ;;
   -sd)
       SUPERFLLISTDIR=${2}
       shift 2
       ;;
   -scb)
       SCIENCEBASEDIR=${2}
       shift 2
       ;;
   -stb)
       STANDARDBASEDIR=${2}
       shift 2
       ;;
   -su)
       if [ "${2}" = "N" ]; then
         SUPERFLAT=""
       fi
       shift 2
       ;;
   -t)
       FLATTYPE=${2}
       shift 2
       ;;
   -tracks)
       TRACKS=1  
       shift 1
       ;;
   -www)
       WEBDIRROOT=${2}
       shift 2
       ;;
    *)
       GoOn=1
       ;;
   esac
done

# sanity checks:
if [ "${MODE}_A" = "_A" ]; then
  echo "nothing to do: exiting!"
  exit 1
fi

if [ "${TRACKS}" -eq 1 ] && [ "${NOTRACKS}" -eq 1 ]; then
  echo "TRACKS and NOTRACKS are both set! Please correct!"
  echo "Exiting!"
  exit 1
fi

export BASEDIR
export FILTER
export RUN
export MD=${BASEDIR}/${SURVEY}/${FILTER}/${RUN}
export TEMPDIR=`pwd`

# number of processors on the machine (assumes Linux):
NPROC=1
if [ -f /proc/cpuinfo ]; then
  NPROC=`grep processor /proc/cpuinfo | wc -l`

  # just limit the number of parallel channels to avoid
  # too heavy I/O:
  if [ ${NPROC} -gt 16 ]; then
    NPROC=16
  fi
fi

# root directory for the sets
if [ "${SETBASEDIR}" = "" ]; then
  SETDIR=${BASEDIR}/${SURVEY}
else
  SETDIR=${SETBASEDIR}/${SURVEY}
fi

REGDIR=${BASEDIR}/regs/
REGPRESENT=0  # if '1' then region files for the run under 
              # consideration are present. The variable is finally
              # set in the TESTREG MODE below.

# root directory for WEB pages:
WEBDIR=${WEBDIRROOT}/${SURVEY}WWW

# if fringing is done the endings of the
# science fromes in the individual set directories
# are OFC(S)F instead of OFC(S):
if [ "${FILTER}" = "i_SDSS" ] ||\
   [ "${FILTER}" = "z_SDSS" ]; then
  FRINGING="F"
else
  FRINGING=""
fi

export MAGMIN=-100  # dummy for the automatic absolute photometric
                    # calibration mode
export MAGMAX=100   #      "

PHOTCAT=`pwd`/STRIPE82.cat
#PHOTCAT=/export/euclid2_1/terben/reduce_KIDS/SLOAN_KIDS_standards.cat
# default values for color index, extinction coefficient
# and color term (absolute photometry) in the different
# WFI filters
#
if [ "${FILTER}" = "g_SDSS" ]; then
  export COLOR=gmr
  export EXT=-0.11
  export COLCOEFF=0.1
  export MAGMIN=24.6  # for the automatic absolute photometric
  export MAGMAX=25.0
fi

if [ "${FILTER}" = "i_SDSS" ]; then
  export COLOR=rmi
  export EXT=-0.04
  export COLCOEFF=0.00
  export MAGMIN=24.0  # for the automatic absolute photometric
  export MAGMAX=25.0
fi

if [ "${FILTER}" = "r_SDSS" ]; then
  export COLOR=gmr
  export EXT=-0.1
  export COLCOEFF=0.05
  export MAGMIN=23.0  # for the automatic absolute photometric
  export MAGMAX=26.0
fi

if [ "${FILTER}" = "u_SDSS" ]; then
  export COLOR=umg
  export EXT=-0.4
  export COLCOEFF=-0.02
  export MAGMIN=23.0  # for the automatic absolute photometric
  export MAGMAX=26.0
fi

#
# default values for thresshholds of global weights:
#
export FLATMIN=0.7
export FLATMAX=1.3
export SCIENCEMIN=0.85
export SCIENCEMAX=1.15
export DARKMIN=-7
export DARKMAX=10

#
# default values for threshholds of science images
# to be marked as bad:
export SCIENCEBADMIN=120
export SCIENCEBADMAX=5000

export SUPERFLATTHRESH=50

#
# modification of thressholds for specific filters:
# OMEGACAM has extremly low counts in u and we need to adapt
# some linits!!
if [ "${FILTER}" = "u_SDSS" ]; then
  SUPERFLATTHRESH=3
  SCIENCEMIN=0.5
  SCIENCEMAX=1.5
fi

#
## Here starts the real work !!!!
#

# preparatory work (splitting of images; checking of integrity);
# note that BIAS and DARKS are typically processed by the script
# 'doall_bias_dark_OMEGACAM_single.sh', and FLATS by
# 'doall_flat_OMEGACAM_single.sh'

# unzip RICE compressed data:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "UNZIP" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 6)}'`
    ./uncompress.sh -s ${MD}/${IMTYPE}_${FILTER}/ORIGINALS \
                    -p ${NPROC} -m UNCOMPRESS
  fi
done

# split data:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "PREPARE" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 8)}'`  
    if [ -d ${MD}/${IMTYPE}_${FILTER} ]; then
      ./process_split_OMEGACAM_eclipse.sh -md ${MD} -sd ${IMTYPE}_${FILTER}

      LIMITS="${SCIENCEBADMIN} ${SCIENCEBADMAX}"
      if [ "${IMTYPE}" != "SCIENCE" ]; then
        # limits for valid standard or short exposed science
        # observations:  
        LIMITS="50 3500"  
      fi

      ./check_files.sh ${MD} ${IMTYPE}_${FILTER} "" ${LIMITS}

      test -d ${MD}/${IMTYPE}_${FILTER}/CT && \
       rm -rf ${MD}/${IMTYPE}_${FILTER}/CT

##       # correct science data for crosstalk:
##       if [[ "${IMTYPE}" =~ "SCIENCE" ]]; then
##         # if the coefficients have already been determined at some
##         # earlier stage we do not need to repeat this:
##         if [ ! -s ${CTLISTDIR}/ct_coeffs_${IMTYPE}_${FILTER}_${RUN} ]; then
##           ./create_crosstalk_coefficients.sh ${MD} ${IMTYPE}_${FILTER} \
##               OBSTART 0.01 ${CTLISTDIR}/ct_coeffs_${IMTYPE}_${FILTER}_${RUN}
##         fi
##         ./apply_crosstalk.sh ${MD} ${IMTYPE}_${FILTER} \
##               ${CTLISTDIR}/ct_coeffs_${IMTYPE}_${FILTER}_${RUN} \
##               ${MD}/${IMTYPE}_${FILTER}/CT ""
##       fi
##       if [ -d ${MD}/${IMTYPE}_${FILTER}/CT ]; then
##         mv ${MD}/${IMTYPE}_${FILTER}/CT/*fits ${MD}/${IMTYPE}_${FILTER}
##         rmdir ${MD}/${IMTYPE}_${FILTER}/CT
##       fi

      # by default remove the BADMODE directory which contains
      # images that we do not further process:
      test -d ${MD}/${IMTYPE}_${FILTER}/BADMODE && \
       rm -rf ${MD}/${IMTYPE}_${FILTER}/BADMODE
    fi
  fi
done

# link/copy processed BIAS frames:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "BIAS" ]]; then
    COPYMODE=`echo ${mode} | awk '{print substr($0, 5)}'`
    if [ "${BIASBASEDIR}" != "" ]; then
      test -d ${MD}/BIAS || mkdir ${MD}/BIAS
      if [ "${COPYMODE}" = "LINK" ]; then
        ln -s ${BIASBASEDIR}/${RUN}/BIAS/BIAS_*fits ${MD}/BIAS
      else
        cp ${BIASBASEDIR}/${RUN}/BIAS/BIAS_*fits ${MD}/BIAS
      fi
    fi
  fi
done

# link/copy processed DARK frames:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "DARK" ]]; then
    COPYMODE=`echo ${mode} | awk '{print substr($0, 5)}'`
    if [ "${DARKBASEDIR}" != "" ]; then
      test -d ${MD}/DARK || mkdir ${MD}/DARK
      if [ "${COPYMODE}" = "LINK" ]; then
        ln -s ${DARKBASEDIR}/${RUN}/DARK/DARK_*fits ${MD}/DARK
      else
        cp ${DARKBASEDIR}/${RUN}/DARK/DARK_*fits ${MD}/DARK
      fi
    fi
  fi
done

# link/copy processed FLAT frames:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "FLAT" ]]; then
    COPYMODE=`echo ${mode} | awk '{print substr($0, 5)}'`
    if [ "${FLATBASEDIR}" != "" ]; then
      test -d ${MD}/${FLATTYPE}_${FILTER} || mkdir ${MD}/${FLATTYPE}_${FILTER}
      if [ "${COPYMODE}" = "LINK" ]; then
        ln -s ${FLATBASEDIR}/${RUN}/${FLATTYPE}_${FILTER}/${FLATTYPE}_${FILTER}_[0-9]*fits \
        ${MD}/${FLATTYPE}_${FILTER}
      else
        cp ${FLATBASEDIR}/${RUN}/${FLATTYPE}_${FILTER}/${FLATTYPE}_${FILTER}_[0-9]*fits \
        ${MD}/${FLATTYPE}_${FILTER}
      fi
    fi
  fi
done

# link/copy processed SCIENCE master frames:
for mode in ${MODE}
do
  # we cannot use 'if [[ "${mode}" =~ "SCIENCE" ]]; then' in the following
  # because of other existing SCIENCE modes:
  if [ "${mode}" = "SCIENCELINK" ] || [ "${mode}" = "SCIENCECOPY" ]; then
    COPYMODE=`echo ${mode} | awk '{print substr($0, 8)}'`
    if [ "${SCIENCEBASEDIR}" != "" ]; then
      test -d ${MD}/SCIENCE_${FILTER} || mkdir ${MD}/SCIENCE_${FILTER}
      if [ "${COPYMODE}" = "LINK" ]; then
        ln -s ${SCIENCEBASEDIR}/${RUN}/SCIENCE_${FILTER}/SCIENCE_${FILTER}_[1-9].fits \
        ${MD}/SCIENCE_${FILTER}
        ln -s ${SCIENCEBASEDIR}/${RUN}/SCIENCE_${FILTER}/SCIENCE_${FILTER}_[1-9][0-9].fits \
        ${MD}/SCIENCE_${FILTER}
      else
        cp ${SCIENCEBASEDIR}/${RUN}/SCIENCE_${FILTER}/SCIENCE_${FILTER}_[1-9].fits \
        ${MD}/SCIENCE_${FILTER}
        cp ${SCIENCEBASEDIR}/${RUN}/SCIENCE_${FILTER}/SCIENCE_${FILTER}_[1-9][0-9].fits \
        ${MD}/SCIENCE_${FILTER}
      fi
    fi
  fi
done

# copy standardstar photometric info:
for mode in ${MODE}
do
  if [ "${mode}" = "PHOTINFOCOPY" ]; then
    if [ "${STANDARDBASEDIR}" != "" ]; then
      test -d ${MD}/STANDARD_${FILTER}/calib ||\
        mkdir -p ${MD}/STANDARD_${FILTER}/calib
      cp -r ${STANDARDBASEDIR}/${RUN}/STANDARD_${FILTER}/calib \
        ${MD}/STANDARD_${FILTER}
    fi
  fi
done

# preprocessing of science frames
for mode in ${MODE}
do
  if [ "${mode}" = "SCIENCE" ]; then
    # header update
    ./parallel_manager.sh ./create_runid_para.sh ${MD} \
        SCIENCE_${FILTER} "." ${RUN}

    # preprocessing

    # check for superflat exclusion lists:
    SUPERFLLIST=0
    if [ -f ${SUPERFLLISTDIR}/superflat_exclusion_${FILTER}_${RUN} ]
    then
      cp ${SUPERFLLISTDIR}/superflat_exclusion_${FILTER}_${RUN} \
         ${MD}/superflat_exclusion
      SUPERFLLIST=1
    fi

    # initial science processing:
    if [ "${FRINGING}" = "F" ]; then
      ./parallel_manager.sh ./process_science_eclipse_para.sh ${MD} BIAS \
          ${FLATTYPE}_${FILTER} \
          SCIENCE_${FILTER} NORESCALE FRINGE SUPERTEST
    else
      if [ "${FILTER}" = "u_SDSS" ] || [ "${SUPERFLAT}_A" = "_A" ]; then
        ./parallel_manager.sh ./process_science_eclipse_para.sh ${MD} \
          BIAS ${FLATTYPE}_${FILTER} \
          SCIENCE_${FILTER} RESCALE NOFRINGE SUPERTEST
      else
        ./parallel_manager.sh ./process_science_eclipse_para.sh ${MD} \
          BIAS ${FLATTYPE}_${FILTER} \
          SCIENCE_${FILTER} NORESCALE NOFRINGE SUPERTEST
      fi
    fi

    # by default remove splitted images: 
    if [ ${CLEANUP} -eq 1 ]; then
      test -d ${MD}/SCIENCE_${FILTER}/SPLIT_IMAGES && \
       rm -rf ${MD}/SCIENCE_${FILTER}/SPLIT_IMAGES
    fi
    
    # copy superflatlists if necessary:
    if [ ${SUPERFLLIST} -eq 0 ]; then
      TEMP=""
      for i in ${MD}/superflat_exclusion_*
      do
        TEMP="${TEMP} ${i}"
      done

      if [ "${TEMP}" != "" ]; then
        test -d ${SUPERFLLISTDIR} || mkdir -p ${SUPERFLLISTDIR}
  
        cat ${TEMP} > \
          ${SUPERFLLISTDIR}/superflat_exclusion_${FILTER}_${RUN}
      fi
    fi

    rm ${MD}/superflat_exclusion_*

    ./parallel_manager.sh ./create_illumfringe_para.sh ${MD} SCIENCE_${FILTER}
    #./resolvelinks.sh ${MD}/SCIENCE_${FILTER} SCIENCE_${FILTER} illum

    # do superflat correction:
    if [ "${SUPERFLAT}" = "S" ]; then
      if [ "${FILTER}" = "u_SDSS" ]; then
         ./parallel_manager.sh ./process_science_illum_eclipse_para.sh ${MD}\
                   SCIENCE_${FILTER} NORESCALE ILLUM \
                   ${SUPERFLATTHRESH}
      else
         ./parallel_manager.sh ./process_science_illum_eclipse_para.sh ${MD}\
                   SCIENCE_${FILTER} RESCALE ILLUM \
                   ${SUPERFLATTHRESH}
      fi

      # by default remove OFC images and sky-subtracted versions
      # that are only used for superlfat creation: 
      if [ ${CLEANUP} -eq 1 ]; then
        test -d ${MD}/SCIENCE_${FILTER}/OFC_IMAGES && \
         rm -rf ${MD}/SCIENCE_${FILTER}/OFC_IMAGES
        test -d ${MD}/SCIENCE_${FILTER}/SUB_IMAGES && \
         rm -rf ${MD}/SCIENCE_${FILTER}/SUB_IMAGES
      fi
    fi

    # do fringe correction:
    if [ "${FRINGING}" = "F" ]; then
      ./parallel_manager.sh ./process_science_fringe_eclipse_para.sh ${MD} \
                              SCIENCE_${FILTER}

      # by default remove OFC/OFCS images: 
      if [ ${CLEANUP} -eq 1 ]; then
        test -d ${MD}/SCIENCE_${FILTER}/OFC_IMAGES && \
         rm -rf ${MD}/SCIENCE_${FILTER}/OFC_IMAGES
        test -d ${MD}/SCIENCE_${FILTER}/OFCS_IMAGES && \
         rm -rf ${MD}/SCIENCE_${FILTER}/OFCS_IMAGES
      fi
    fi
  fi
done


# Preprocess standards and other, short exposed science frames
for mode in ${MODE}
do
  if [ "${mode}" = "STANDARD" ] || [ "${mode}" = "SCIENCESHORT" ]; then
    # Preprocess the standards without creating superflat again
    if [ -d ${MD}/${mode}_${FILTER} ]; then
      if [ "${FILTER}" = "u_SDSS" ]; then
        ./parallel_manager.sh ./process_standard_eclipse_para.sh ${MD} BIAS \
            ${FLATTYPE}_${FILTER} SCIENCE_${FILTER} ${mode}_${FILTER} RESCALE
      else
        ./parallel_manager.sh ./process_standard_eclipse_para.sh ${MD} BIAS \
            ${FLATTYPE}_${FILTER} SCIENCE_${FILTER} \
            ${mode}_${FILTER} NORESCALE
      fi

      # by default remove splitted images: 
      if [ ${CLEANUP} -eq 1 ]; then
        test -d ${MD}/${mode}_${FILTER}/SPLIT_IMAGES && \
         rm -rf ${MD}/${mode}_${FILTER}/SPLIT_IMAGES
      fi

      ./parallel_manager.sh ./create_illumfringe_para.sh ${MD} \
          ${mode}_${FILTER}

      if [ "${FILTER}" = "u_SDSS" ]; then
        ./parallel_manager.sh ./process_science_illum_eclipse_para.sh ${MD} \
            ${mode}_${FILTER} NORESCALE ILLUM ${SUPERFLATTHRESH}

      else
        ./parallel_manager.sh ./process_science_illum_eclipse_para.sh ${MD} \
            ${mode}_${FILTER} RESCALE ILLUM ${SUPERFLATTHRESH}

      fi

      # by default remove OFC images: 
      if [ ${CLEANUP} -eq 1 ]; then
        test -d ${MD}/${mode}_${FILTER}/OFC_IMAGES && \
         rm -rf ${MD}/${mode}_${FILTER}/OFC_IMAGES
      fi

      if [ "${FRINGING}" = "F" ]; then
        ./parallel_manager.sh ./process_science_fringe_eclipse_para.sh ${MD}\
            ${mode}_${FILTER}

        # by default remove OFCS images: 
        if [ ${CLEANUP} -eq 1 ]; then
          test -d ${MD}/${mode}_${FILTER}/OFCS_IMAGES && \
           rm -rf ${MD}/${mode}_${FILTER}/OFCS_IMAGES
        fi
      fi
    fi
  fi
done

# create binned mosaics
for mode in ${MODE}
do
  if [[ "${mode}" =~ "MOSAICS" ]]; then
     IMTYPE=`echo ${mode} | awk '{print substr($0, 8)}'`
     ./create_binnedmosaics_exposurepara.sh ${MD} ${IMTYPE}_${FILTER} "*" \
         OFC${SUPERFLAT}${FRINGING} 8 -32
  fi
done

# create global weight images
for mode in ${MODE}
do
  if [ "${mode}" = "GLOBALWEIGHTS" ]; then
    # delete an old WEIGHTS directory:
    test -d ${MD}/WEIGHTS && rm -rf ${MD}/WEIGHTS
    test -d ${MD}/SCIENCE_${FILTER}_norm && \
     rm -rf ${MD}/SCIENCE_${FILTER}_norm
    test -d ${MD}/${FLATTYPE}_${FILTER}_norm && \
     rm -rf ${MD}/${FLATTYPE}_${FILTER}_norm
  
    ./parallel_manager.sh create_norm_para.sh ${MD} ${FLATTYPE}_${FILTER}
    ./parallel_manager.sh create_norm_para.sh ${MD} SCIENCE_${FILTER}

    ./parallel_manager.sh ./create_global_weights_flags_para.sh ${MD} \
      NOLINK ${FLATTYPE}_${FILTER}_norm ${FLATMIN} ${FLATMAX} DARK \
      ${DARKMIN} ${DARKMAX} SCIENCE_${FILTER}_norm ${SCIENCEMIN} ${SCIENCEMAX}
  fi
done

# check for the existence of region files
for mode in ${MODE}
do
  if [ "${mode}" = "TESTREG" ]; then
    if [ ${TRACKS} -eq 0 ]; then
      if [ ${NOTRACKS} -eq 1 ]; then
        REGPRESENT=1
      else
        # REG files are 'defined' to not be present if:
        # - The corresponding reg file directory does not exist
        # - The corresponding reg file directory does exist but
        #   is empty
        # - There are only reg files for short exposures (less than
        #   100 sec.)
        REGPRESENT=0
        if [ -d ${REGDIR}/${FILTER}/${RUN} ]; then
          for REGFILE in `find ${REGDIR}/${FILTER}/${RUN} -name \*reg`
          do
            # check exptime of all FITS files of which reg files
            # exist
            FITSFILE=${MD}/SCIENCE_${FILTER}/`basename ${REGFILE} .reg`\
                     OFC${SUPERFLAT}${FRINGING}.fits

            if [ -f ${FITSFILE} ]; then
              SCIENCE_FILE=`dfits ${FITSFILE} | fitsort -d EXPTIME\
                            | awk '{if($2 > 100.0) { print 1} else {print 0}}'`

              if [ ${SCIENCE_FILE} -eq 1 ]; then
                REGPRESENT=1
              fi
            fi
          done
        fi
      fi
    fi
  fi
done

# create weight images (applicable to SCIENCE and SCIENCESHORT exposures)
for mode in ${MODE}
do
  if [[ "${mode}" =~ "WEIGHTS" ]] && [ "${mode}" != "GLOBALWEIGHTS" ]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 8)}'`

    if [ -z ${IMTYPE} ]; then
      echo "Unknown image type in WEIGHTS mode! Exiting!"
      exit 1;  
    fi

    # use region files if present:
    if [ ${REGPRESENT} -eq 1 ]; then
      test -d ${MD}/${IMTYPE}_${FILTER}/reg ||\
           mkdir ${MD}/${IMTYPE}_${FILTER}/reg
      cp ${REGDIR}/${FILTER}/${RUN}/*reg ${MD}/${IMTYPE}_${FILTER}/reg
    fi
    ./parallel_manager.sh ./create_weights_flags_para.sh ${MD} \
            ${IMTYPE}_${FILTER} OFC${SUPERFLAT}${FRINGING} WEIGHTS_FLAGS
  fi
done

# create Satellite TRACK images
for mode in ${MODE}
do
  if [ "${mode}" = "STRACK" ]; then
    # do nothing if region files are already available (REGPRESENT=1
    # in this case):
    if [ ${REGPRESENT} -eq 0 ]; then
      # because the satellite track detection creates very large files
      # we split up the process in junks of 'NFILES' files.
      NSTRACKFILES=5

      find ${MD}/SCIENCE_${FILTER} -name \*OFC${SUPERFLAT}${FRINGING}.fits |\
        awk -F/ '{print $NF}' | awk -F_ '{print $1}' | sort | uniq |\
        awk 'BEGIN { n = 1; name = "'${TEMPDIR}'/" "SEG_images_" n "_'$$'"; } {
               print $0 >> name;
               if(NR == n * '${NSTRACKFILES}') {
                 n++;
                 name = "SEG_images_" n "_'$$'";
               }
             }'

      for FILE in ${TEMPDIR}/SEG_images*_$$
      do
        # Create the hough images
        ./parallel_manager.sh ./create_hough_para.sh \
          ${MD} SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING} ${FILE}

        # Detect the tracks
        ./create_strack_mask_exposurepara.sh \
          ${MD} SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING} ${FILE}

        # delete not needes stuff from the track detection.
        # Especially the 'SEG' images are very large:
        ./cleanfiles.sh ${MD}/SCIENCE_${FILTER} '*' \
           "OFC${SUPERFLAT}${FRINGING}_SEG."
        ./cleanfiles.sh ${MD}/SCIENCE_${FILTER} '*' \
           "OFC${SUPERFLAT}${FRINGING}_SEG_hSN."
        ./cleanfiles.sh ${MD}/SCIENCE_${FILTER} '*' \
           "OFC${SUPERFLAT}${FRINGING}_SEG.fits_weight_ref."
        find ${MD}/SCIENCE_${FILTER} -maxdepth 1 -name default\* | xargs rm
        find ${MD}/SCIENCE_${FILTER} -maxdepth 1 -name tmp\* | xargs rm
        find ${MD}/SCIENCE_${FILTER}/cat -maxdepth 1 \
             -name seeing\*cat | xargs rm
        find ${MD}/SCIENCE_${FILTER} -maxdepth 1 -name \*LOG | xargs rm
        find ${MD}/SCIENCE_${FILTER} -maxdepth 1 -name log\* | xargs rm
        find ${MD}/SCIENCE_${FILTER} -maxdepth 1 -name HTweight\*fits\* \
             | xargs rm

        # copy region files to appropriate subdirs of REGDIR
        # and to a reg subdirectory:
        test -d ${MD}/SCIENCE_${FILTER}/reg || \
          mkdir ${MD}/SCIENCE_${FILTER}/reg
        test -d ${REGDIR}/${FILTER}/${RUN}  || \
          mkdir -p ${REGDIR}/${FILTER}/${RUN}

        FILES=`find ${MD}/SCIENCE_${FILTER}/STRACK \
                    -name \*OFC${SUPERFLAT}${FRINGING}.reg`
          
        for FILE in ${FILES}
        do
          BASE=`basename ${FILE} OFC${SUPERFLAT}${FRINGING}.reg`
          cp ${FILE} ${MD}/SCIENCE_${FILTER}/reg/${BASE}.reg
          cp ${FILE} ${REGDIR}/${FILTER}/${RUN}/${BASE}.reg
        done
      done
      rm ${TEMPDIR}/SEG_images*_$$
    fi
  fi
done

# create initial weight images
for mode in ${MODE}
do
  if [ "${mode}" = "REGIONWEIGHTS" ]; then
    # we definitely do not need to execute this mode  
    # if the user forces notrcks:
    if [ ${NOTRACKS} -ne 1 ]; then  
      if [ -d  ${MD}/SCIENCE_${FILTER}/reg ]; then
        # first delete old weight files:
        ./cleanfiles.sh ${MD}/WEIGHTS '*' "OFC${SUPERFLAT}${FRINGING}.weight."
        ./cleanfiles.sh ${MD}/WEIGHTS '*' "OFC${SUPERFLAT}${FRINGING}.flag."

        ./parallel_manager.sh ./create_weights_flags_para.sh ${MD} \
              SCIENCE_${FILTER} "OFC${SUPERFLAT}${FRINGING}" WEIGHTS_FLAGS
      fi
    fi
  fi
done

# 'correct' weight/flag images for too agressive stellar masking
# in very good seeing conditions (cosmic ray mask problem)
for mode in ${MODE}
do
  if [ "${mode}" = "CORRECTWEIGHTS" ]; then
    # first create catalogues WITHOUT the weights:
    # as we aim for high S/N stars we omit details on thresholds
    # for different filters;
    #
    # Only do it for g, r, i, y and z:
    if [ "${FILTER}" = "g_SDSS" ] || [ "${FILTER}" = "r_SDSS" ] ||\
       [ "${FILTER}" = "i_SDSS" ] || [ "${FILTER}" = "z_SDSS" ]; then
      ./parallel_manager.sh ./create_astromcats_weights_para.sh ${MD} \
          SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING}\
          cat_maskcorr NONE NONE NONE 10 200 3

      # we 'abuse' the check PSF script to create KSB catalogues:
      ./parallel_manager.sh ./check_science_PSF_para.sh ${MD} \
          SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING}\
          cat_maskcorr

      ./merge_sex_ksb.sh ${MD} SCIENCE_${FILTER} \
          OFC${SUPERFLAT}${FRINGING}\
          cat_maskcorr

      # The '4.1' in the following call is the fixed radius around stars
      # in which weights/flags are 'fixed'.
      ./parallel_manager.sh ./fix_stars_weights_para.sh \
          ${MD} SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING}\
          cat_maskcorr WEIGHTS OFC${SUPERFLAT}${FRINGING}.weight \
          OFC${SUPERFLAT}${FRINGING}.flag 4.1 2.0

      # by default remove the old, uncorrected weights: 
      if [ ${CLEANUP} -eq 1 ]; then
        test -d ${MD}/WEIGHTS/orig && \
         rm -rf ${MD}/WEIGHTS/orig
        test -d ${MD}/SCIENCE_${FILTER}/cat_maskcorr && \
         rm -rf ${MD}/SCIENCE_${FILTER}/cat_maskcorr
      fi
    fi
  fi
done

# correct the weights of chip 25 for crosstalk effects by chip 26:
# for the moment we only do this for 'science' KIDS images:
for mode in ${MODE}
do
  if [ "${mode}" = "CROSSTALKWEIGHTS" ]; then
    ./create_OMEGACAM_ct_weights.sh ${MD} SCIENCE_${FILTER} \
       OFC${SUPERFLAT}${FRINGING} WEIGHTS \
       "25 26" "26 25" "25 27" "27 25" "26 27" "27 26"
  fi
done

# create weight images, astrometrically calibrate standard
# images. Do photometric calibration
for mode in ${MODE}
do
  if [ "${mode}" = "ABSPHOTOM" ]; then
    if [ -d ${MD}/STANDARD_${FILTER} ]; then
      ./parallel_manager.sh ./create_weights_flags_para.sh ${MD} \
         STANDARD_${FILTER} OFC${SUPERFLAT}${FRINGING} WEIGHTS_FLAGS

      ./parallel_manager.sh ./create_astromcats_weights_para.sh \
         ${MD} STANDARD_${FILTER} OFC${SUPERFLAT}${FRINGING} cat \
         WEIGHTS OFC${SUPERFLAT}${FRINGING}.weight \
         OFC${SUPERFLAT}${FRINGING}.flag 5

      ./create_scamp_astrom_photom.sh ${MD} STANDARD_${FILTER}\
         OFC${SUPERFLAT}${FRINGING} 2.0 2MASS

      ./parallel_manager.sh ./create_stdphotom_prepare_para.sh ${MD} \
         STANDARD_${FILTER} OFC${SUPERFLAT}${FRINGING} 2MASS ${PHOTCAT}

      # filter name in PHOTCAT. It is just the first letter of
      # the complete FILTER name:
      FILTNAMSTAND=`echo ${FILTER} | awk '{print substr($0, 1, 1)}'`
      ./create_abs_photo_info.sh ${MD} STANDARD_${FILTER} \
         SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING} \
         ${FILTER} ${FILTNAMSTAND} ${COLOR} ${EXT} \
         ${COLCOEFF} RUNCALIB\
         AUTOMATIC ${MAGMIN} ${MAGMAX}
    fi
    ./create_zp_correct_header.sh ${MD} SCIENCE_${FILTER} \
         OFC${SUPERFLAT}${FRINGING}
  fi
done

# transfer photometric solution from the main science frames also to the
# shoer science observations:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "PHOTOM" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 7)}'`

    if [ -f ${MD}/STANDARD_${FILTER}/calib/night_0_${FILTER}_result.asc ]
    then
      ./create_zp_correct_header.sh ${MD} ${IMTYPE}_${FILTER} \
         OFC${SUPERFLAT}${FRINGING} \
         ${MD}/STANDARD_${FILTER}/calib/night_0_${FILTER}_result.asc 2 RUN
    else
      echo "File ${MD}/STANDARD_${FILTER}/calib/night_0_${FILTER}_result.asc missing."  
      echo "We mark images as non-photometrically calibrated!"
      ./create_zp_correct_header.sh ${MD} ${IMTYPE}_${FILTER} \
         OFC${SUPERFLAT}${FRINGING}
    fi
  fi
done

# create astrometric catalogues (this is mainly to be able to create
# WEB pages on the run basis)
for mode in ${MODE}
do
  if [[ "${mode}" =~ "SINGLEASTROM" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 13)}'`

    # delete old catalogues if present:
    test -d ${MD}/${IMTYPE}_${FILTER}/cat && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/cat

    # default threshhold is 5 pixels with 5 sigma (for 'normal' data)
    THRESH=5

    if [ "${FILTER}" = "u_SDSS" ]; then
      THRESH=4
    fi

    ./parallel_manager.sh ./create_astromcats_weights_para.sh ${MD} \
           ${IMTYPE}_${FILTER} OFC${SUPERFLAT}${FRINGING} cat\
           WEIGHTS OFC${SUPERFLAT}${FRINGING}.weight NONE ${THRESH}

    #
    # to make things simple for the WEB page creation we create
    # PSF checkplots for all colours
    ./parallel_manager.sh ./check_science_PSF_para.sh ${MD} \
            ${IMTYPE}_${FILTER} OFC${SUPERFLAT}${FRINGING} cat

    ./check_science_PSF_plot.sh ${MD} ${IMTYPE}_${FILTER} \
        OFC${SUPERFLAT}${FRINGING} cat

    ./merge_sex_ksb.sh ${MD} ${IMTYPE}_${FILTER} \
        OFC${SUPERFLAT}${FRINGING} cat

    ./create_stats_table.sh ${MD} ${IMTYPE}_${FILTER} \
        OFC${SUPERFLAT}${FRINGING}
  fi
done

# Create web pages for SCIENCE frames
for mode in ${MODE}
do
  if [[ "${mode}" =~ "RUNPAGES" ]]; then
    WWWTYPE=`echo ${mode} | awk '{print substr($0, 9)}'`
      
    if [ "${WWWTYPE}" = "SCIENCE" ]; then
      WWWMODE=runscience  
    else
      WWWMODE=runscienceshort
    fi

    ./pipelog.py -a ${WWWMODE} -f ${FILTER} -o ${WEBDIR} -p "RUN"\
            -b ${BASEDIR}/${SURVEY} -r ${RUN} -e OFC${SUPERFLAT}${FRINGING}
    ./pipelog.py -a runindex -f ${FILTER} -o ${WEBDIR} -p "RUN"\
            -b ${BASEDIR}/${SURVEY} -r ${RUN}
  fi
done

# subtract the sky from science images:
for mode in ${MODE}
do
  if [ "${mode}" = "SKYSUB" ]; then
    # add sky backgroundvalue to the image headers - we need this in the
    # single images for lensfit:
    ./parallel_manager.sh ./create_addbackgrkey_para.sh ${MD}\
          SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING}

    ./parallel_manager.sh ./create_skysub_para.sh  ${MD} SCIENCE_${FILTER}\
          OFC${SUPERFLAT}${FRINGING} ".sub" TWOPASS_X_COLLAPSE
  fi
done


# Update SCIENCE WWW pages with ZPs
for mode in ${MODE}
do
  if [ "${mode}" = "UPDATERUNSCIENCEPAGES" ]; then
    ./pipelog.py -a runscience -f ${FILTER} -o ${WEBDIR} -b ${BASEDIR} \
        -r ${RUN} -e OFC${SUPERFLAT}${FRINGING} -p "RUN" -u
  fi
done

# distribute current run
for mode in ${MODE}
do
  if [ "${mode}" = "RUNDISTRIBUTE" ]; then
    ./parallel_manager.sh ./link_sets_OMEGACAM_para.sh ${MD} \
        SCIENCE_${FILTER} OFC${SUPERFLAT}${FRINGING} ${FILTER} \
        ${SETDIR} ${SURVEY}
  fi
done

# clean files from the SCIENCE etc. directories to bring them
# to the original state:
for mode in ${MODE}
do
  if [[ "${mode}" =~ "CLEAN" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 6)}'`
    if [ "${FRINGING}" = "F" ]; then
      ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/OFCS_IMAGES OMEGA OFCS.
      test -d ${MD}/${IMTYPE}_${FILTER}/OFCS_IMAGES && \
        rmdir ${MD}/${IMTYPE}_${FILTER}/OFCS_IMAGES

      ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/ OMEGA OFCSF.
    else
      ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/ OMEGA OFCS.
    fi
    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/ "" fringe.
    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/ "" illum.
    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/ "" .

    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/OFC_IMAGES OMEGA OFC.
      test -d ${MD}/${IMTYPE}_${FILTER}/OFC_IMAGES && \
        rmdir ${MD}/${IMTYPE}_${FILTER}/OFC_IMAGES
    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/SUB_IMAGES OMEGA "OFC_sub."
      test -d ${MD}/${IMTYPE}_${FILTER}/SUB_IMAGES && \
        rmdir ${MD}/${IMTYPE}_${FILTER}/SUB_IMAGES
    ./cleanfiles.sh ${MD}/${IMTYPE}_${FILTER}/SPLIT_IMAGES OMEGA "."
      test -d ${MD}/${IMTYPE}_${FILTER}/SPLIT_IMAGES && \
        rmdir ${MD}/${IMTYPE}_${FILTER}/SPLIT_IMAGES

    test -d ${MD}/${IMTYPE}_${FILTER}/BADMODE && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/BADMODE
    test -d ${MD}/${IMTYPE}_${FILTER}/cat && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/cat
#    test -d ${MD}/${IMTYPE}_${FILTER}/calib && \
#     rm -rf ${MD}/${IMTYPE}_${FILTER}/calib
    test -d ${MD}/${IMTYPE}_${FILTER}/cat_maskcorr && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/cat_maskcorr
    test -d ${MD}/${IMTYPE}_${FILTER}/reg && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/reg
    test -d ${MD}/${IMTYPE}_${FILTER}/STRACK && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/STRACK
    test -d ${MD}/${IMTYPE}_${FILTER}/BINNED && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/BINNED
    test -d ${MD}/${IMTYPE}_${FILTER}/headers_scamp* && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/headers_scamp*
    test -d ${MD}/${IMTYPE}_${FILTER}/astrom_photom_scamp* && \
     rm -rf ${MD}/${IMTYPE}_${FILTER}/astrom_photom_scamp*
  fi
done

# delete original FITS files; for KIDS we have limks to the .fz
# versions:
for mode in ${MODE}
do
  if [ "${mode}" = "DELETEORIGS" ]; then
    find ${MD}/SCIENCE_${FILTER}/ORIGINALS -name \*fits |\
      xargs -r rm
    find ${MD}/STANDARD_${FILTER}/ORIGINALS -name \*fits |\
      xargs -r rm
    find ${MD}/SCIENCESHORT_${FILTER}/ORIGINALS -name \*fits |\
      xargs -r rm
  fi
done

# delete directories created during RUN processing:
for mode in ${MODE}
do
  if [ "${mode}" = "DELETEDIRS" ]; then
    test -d ${MD}/WEIGHTS && rm -rf ${MD}/WEIGHTS
    test -d ${MD}/SCIENCE_${FILTER}_norm && \
     rm -rf ${MD}/SCIENCE_${FILTER}_norm
    test -d ${MD}/SKYFLAT_${FILTER}_norm && \
     rm -rf ${MD}/SKYFLAT_${FILTER}_norm
    test -d ${MD}/DOMEFLAT_${FILTER}_norm && \
     rm -rf ${MD}/DOMEFLAT_${FILTER}_norm
  fi
done

# illumination correction (calculate illum correction: ILLUMCORRECTION , apply e.g.ILLUMSCIENCE) :
for mode in ${MODE}
do
  if [[ "${mode}" =~ "ILLUM" ]]; then
    IMTYPE=`echo ${mode} | awk '{print substr($0, 6)}'`
    echo $IMTYPE
    FILTNAMSTAND=`echo ${FILTER} | awk '{print substr($0, 1, 1)}'`
    echo $FILTER
    echo $FILTNAMSTAND
    if [ "${IMTYPE}" == "CORRECTION" ]; then
	./illum_correction.sh ${MD} STANDARD_${FILTER} ${FILTER} \
		2 OFC${SUPERFLAT}${FRINGING} ${FILTNAMSTAND} RUNCALIB
    else
	./illum_apply.sh ${MD} ${IMTYPE}_${FILTER} STANDARD_${FILTER} \
		${FILTER} OFC${SUPERFLAT}${FRINGING} RUNCALIB ${NPROC}
    fi
  fi
done
