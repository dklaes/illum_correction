#!/usr/bin/python
# -*- coding: utf-8 -*-
#numpy-1.6.2 required!

#Importing packages
from __future__ import print_function
import numpy as np
import os
import sys
import getopt
import ldac


def filter_elements(data, key, value, condition):
  if (condition == '='):
    mask = (data[key]==value)
  elif (condition == '!='):
    mask = (data[key]!=value)
  elif (condition == '>'):
    mask = (data[key]>value)
  elif (condition == '<'):
    mask = (data[key]<value)
  elif (condition == '>='):
    mask = (data[key]>=value)
  elif (condition == '<='):
    mask = (data[key]<=value)
  else:
    print("Condition not recognized.")
    sys.exit(1)

  return data.filter(mask)
  


def unique_elements(infile, table, key):
  data = ldac.LDACCat(infile)[table]
  unique = np.unique(data[key])
  print(str(unique)[1:-1])

def calcs_before_fitting(infile, outfile, table, external, replace=False):
  data = ldac.LDACCat(infile)[table]
    
  # Getting all information for needed variables.
  # - Mag:     
  # - ZP:
  # - EXT:
  # - AIRMASS:
  # - COLCOEFF:
  # - COLOR:
  Mag = data['Mag']
  MagErr = data['MagErr']
  
  ZP = float(external[0])
  ZPERR = float(external[1])
  
  EXT = float(external[2])
  EXTERR = float(external[3])
  
  AIRMASS = data['AIRMASS']
  AIRMASSERR = 0.0 # not available
  
  COLCOEFF = float(external[4])
  COLCOEFFERR = float(external[5])
  
  colorname = external[6]
  COLOR = data[colorname]
  COLOR1ERR = data[colorname[0] + "mag_err"]
  COLOR2ERR = data[colorname[2] + "mag_err"]
  data[colorname + '_err'] = np.sqrt((COLOR1ERR)**2 + (COLOR2ERR)**2)
  COLORERR = data[colorname + '_err']
  
  filtername = external[7]
  reference = data[filtername]
  reference_err = data[filtername + '_err']
  Xpos_global = data['Xpos_global']
  Ypos_global = data['Ypos_global']
  
  data['MagZP'] = Mag + ZP + EXT*AIRMASS + COLCOEFF*COLOR
  MagZP = data['MagZP']
  
  data['MagZPErr'] = np.sqrt((MagErr)**2 + (ZPERR)**2 + (AIRMASS*EXTERR)**2 + (EXT*AIRMASSERR)**2 + (COLOR*COLCOEFFERR)**2 + (COLCOEFF*COLORERR)**2)
  MagZPErr = data['MagZPErr']
  
  data['Residual'] = MagZP - reference
  Residual = data['Residual']
  
  data['Residual_Err'] = np.sqrt((MagZPErr)**2 + (-reference_err)**2)
  Residual_Err = data['Residual_Err']
  
  data['Xpos_mod'] = 2.0*Xpos_global/PIXXMAX
  Xpos_mod = data['Xpos_mod']
  
  data['Ypos_mod'] = 2.0*Ypos_global/PIXYMAX
  Ypos_mod = data['Ypos_mod']

  data.saveas(outfile, clobber=replace)


opts, args = getopt.getopt(sys.argv[1:], "i:o:t:k:a:v:c:e:", ["input=", "output=", "table=", "key=", "action=", "value=", "condition=", "external="])

infile = outfile = table = key = action = test = value = condition = external = None
for o, a in opts:
    if o in ("-i"):
        infile = a
    elif o in ("-o"):
        outfile = a
    elif o in ("-t"):
        table = a
    elif o in ("-k"):
        key = a
    elif o in ("-a"):
        action = a
    elif o in ("-v"):
	value = int(a)
    elif o in ("-c"):
	condition = a
    elif o in ("-e"):
	external = a.split()

global PIXXMAX
global PIXYMAX
#Reading chip geometry from config file
PIXXMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $1}'").readlines())[0]) * int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
PIXYMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $2}'").readlines())[0]) * int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])




if (action == 'FILTER_ELEMENTS'):
  replace=True
  data = ldac.LDACCat(infile)[table]
  data2 = filter_elements(data, key, value, condition)
  data2.saveas(outfile, clobber=replace)
elif (action == 'UNIQUE_ELEMENTS'):
  unique_elements(infile, table, key)
elif (action == 'PASTE_CATALOGS'):
  ldac.pasteCatalogs(infile, outfile=outfile, table=table, replace=True)
elif (action == 'CALCS_BEFORE_FITTING'):
  # This action contains all calculations that have to be done before fitting.
  # The following argument have to within the "external" string in the following order:
  # zeropoint plus error, extinction coefficient plus error, color coefficient plus error,
  # colorname, filtername.
  calcs_before_fitting(infile, outfile, table, external, replace=True)
elif (action == 'FILTER_PERCENT'):
  # This action contains excluding a lower and upper percent part of the catalog.
  # The following argument have to within the "external" string in the following order:
  # lowerpercent, upperpercent
  replace=True
  data = ldac.LDACCat(infile)[table]
  length = len(data)
  datasorted = np.sort(data[key])
  
  lowerpercent = float(external[0])
  upperpercent = float(external[1])
  
  numberlower = int(length * lowerpercent / 100.0)
  numberupper = int(length * upperpercent / 100.0)
  
  lowervalue = datasorted[numberlower]
  uppervalue = datasorted[length-numberupper]
  
  data2 = filter_elements(data, key, lowervalue, '>')
  data3 = filter_elements(data2, key, uppervalue, '<')
  
  data3.saveas(outfile, clobber=replace)
elif (action == 'FILTER_RESIDUAL'):
  # The following argument have to within the "external" string in the following order:
  # lowercutresabs, uppercutresabs
  replace=True
  data = ldac.LDACCat(infile)[table]
  
  lowercutresabs = float(external[0])
  uppercutresabs = float(external[1])
  
  data2 = filter_elements(data, key, lowercutresabs, '>')
  data3 = filter_elements(data2, key, uppercutresabs, '<')
  
  data3.saveas(outfile, clobber=replace)
elif (action == 'FILTER_RESIDUALMEAN'):
  # The following argument have to within the "external" string in the following order:
  # lowercutresmean, uppercutresmean
  replace=True
  data = ldac.LDACCat(infile)[table]
  
  lowercutresmean = float(external[0])
  uppercutresmean = float(external[1])
  
  mean = np.mean(data[key])
  
  data2 = filter_elements(data, key, lowercutresmean+mean, '>')
  data3 = filter_elements(data2, key, uppercutresmean+mean, '<')
  
  data3.saveas(outfile, clobber=replace)