# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

# ----------------------------------------------------------------
# File Name:           illum_correction_fit.py
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    31.05.2013
# Version:		V1.1
# Description:         Fitting polynomial to data
# ----------------------------------------------------------------

# Changes from V1.0 to V1.1
# - more comments

# $1  main dir

#Importing packages
from __future__ import division, print_function
import numpy as np
import os
import sys
import getopt
import ldac
from scipy.optimize import curve_fit


def delta(x,x0):
  delta = -1.0*abs(x-x0)
  for i in range(len(delta)):
    if int(delta[i])==0:
      delta.itemset(i,1)
    else:
      delta.itemset(i,0)
  return delta

def function_poly2d(xyz, A, B, C, D, E, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21, F22, F23, F24, F25, F26, F27, F28, F29, F30, F31, F32):
    '''
Return a Second order polynomial function:
f(x,y) = A * x ** 2 + B * y ** 2 + C * x * y + D * x + E * y + F[chip]
INPUT
xy
Tuple of numpy arrays containing the x and y coordinates
A, B, C, D, E, F[chip]
Polinomial coefficients
'''
    x, y, z = xyz
    return A * x ** 2 + B * y ** 2 + C * x * y + D * x + E * y + delta(1,z)*F1 + delta(2,z)*F2 + delta(3,z)*F3 + delta(4,z)*F4 + delta(5,z)*F5 + delta(6,z)*F6 + delta(7,z)*F7 + delta(8,z)*F8 + delta(9,z)*F9 + delta(10,z)*F10 + delta(11,z)*F11 + delta(12,z)*F12 + delta(13,z)*F13 + delta(14,z)*F14 + delta(15,z)*F15 + delta(16,z)*F16 + delta(17,z)*F17 + delta(18,z)*F18 + delta(19,z)*F19 + delta(20,z)*F20 + delta(21,z)*F21 + delta(22,z)*F22 + delta(23,z)*F23 + delta(24,z)*F24 + delta(25,z)*F25 + delta(26,z)*F26 + delta(27,z)*F27 + delta(28,z)*F28 + delta(29,z)*F29 + delta(30,z)*F30 + delta(31,z)*F31 + delta(32,z)*F32



#def poly2d_curve_output(best_params, cov,k):
def poly2d_curve_output(best_params, cov):
    pars = ['A', 'B', 'C', 'D', 'E', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23', 'F24', 'F25', 'F26', 'F27', 'F28', 'F29', 'F30', 'F31', 'F32']
    #Diagonal terms of the covariance matrix are the variance of the fitted par
    error = np.sqrt(np.diagonal(cov))
    f = open(path + 'coeffs.txt', 'w')
    for name, value, sig in zip(pars, best_params, error):
      f.write("%s = %2.18f +- %2.18f\n" % (name, value, sig))
    f.close()
    f = open(path + 'coveriance_matrix.txt', 'w')
    f.write("%s" %cov)
    f.close()


#Reading command line arguments
opts, args = getopt.getopt(sys.argv[1:], "i:t:p:", ["input=", "table=", "path="])

infile = table = path = None
for o, a in opts:
    if o in ("-i"):
        infile = a.split()
    elif o in ("-t"):
        table = a
    elif o in ("-p"):
        path = a

#Reading chip geometry from config file
NUMCHIPS = int((os.popen("echo ${NCHIPS} | awk '{print $1}'").readlines())[0])

#Getting x,y coordinates, the residual and the chip number...
data = ldac.LDACCat(infile[0])[table]
x = np.array(data['Xpos_mod'], dtype=np.float64)
y = np.array(data['Ypos_mod'], dtype=np.float64)
eps = np.array(data['Residual'], dtype=np.float64)
chip = np.array(data['CHIP'], dtype=np.int)
sigma = np.array(data['Residual_Err'], dtype=np.float64)


#Now fit...
print("Fitting all chips...")
p0 = np.array((5+NUMCHIPS)*[0.0])
best_par, cov_fit = curve_fit(function_poly2d, (x,y,chip), eps, p0, sigma=sigma)

#Print the output...
poly2d_curve_output(best_par, cov_fit)
print("Fitting complete.")
