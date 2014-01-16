#!/usr/bin/python
# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

#Importing packages
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import os
import sys
import getopt
import ldac
import pyfits


def filter_elements(infile, outfile, table, key, value, condition, replace=False):
  data = ldac.LDACCat(infile)[table]
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

  data = data.filter(mask)
  data.saveas(outfile, clobber=replace)


def unique_elements(infile, table, key):
  data = ldac.LDACCat(infile)[table]
  unique = np.unique(data[key])
  print(str(unique)[1:-1])





opts, args = getopt.getopt(sys.argv[1:], "i:o:t:k:a:v:c:", ["input=", "output=", "table=", "key=", "action=", "value=", "condition="])

infile = outfile = table = key = action = test = value = condition = None
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


if (action == 'FILTER_ELEMENTS'):
  filter_elements(infile, outfile, table, key, value, condition, replace=True)

elif (action == 'UNIQUE_ELEMENTS'):
  unique_elements(infile, table, key)
  
elif (action == 'PASTE_CATALOGS'):
  ldac.pasteCatalogs(infile, outfile=outfile, table=table, replace=True)
