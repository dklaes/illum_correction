#!/usr/bin/python
# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

# ----------------------------------------------------------------
# File Name:           illum_correction_contourplot_fitfunction.py
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    31.05.2013
# Version:		V1.1
# Description:         Create FITS-correction files and corresponding plots
# ----------------------------------------------------------------

# Changes from V1.0 to V1.1
# - more comments
# - Checkplot magnitude / residual dependency added
# - Now running also without X-server
# - old "illum_correction_plot_fitted.py" included and restructed to multiprocessing

# $1  main dir
# $2  lower residual limit (always negative)
# $3  upper residual limit (always positive)

#Importing packages
from __future__ import division, print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
import os
import sys
import multiprocessing
import pyfits


def appearance(xlabel, ylabel, title, xlimits, ylimits, grid, camgrid, plt, sub):
  if (grid == 'grid'):
    plt.grid(True)		# Adding a grid
  elif (grid == 'nogrid'):
    plt.grid(False)
  if (sub == 'nosub'):
    plt.xlabel(xlabel)		# x label
    plt.ylabel(ylabel)		# y label
    plt.title(title)		# Adding a title
    if (xlimits != ''):	# Giving xrange, if empty then auto
	plt.xlim(xlimits)
    if (ylimits != ''):	# Giving yrange, if empty then auto
	plt.ylim(ylimits)
    if (bar == 'colorbar'):
	plt.colorbar()
  elif (sub == 'sub'):
    plt.set_xlabel(xlabel)	# x label
    plt.set_ylabel(ylabel)	# y label
    plt.set_title(title)	# Adding a title
    if (xlimits != ''):	# Giving xrange, if empty then auto
	plt.set_xlim(xlimits)
    if (ylimits != ''):	# Giving yrange, if empty then auto
	plt.set_ylim(ylimits)

  if (camgrid == 'camgrid'):
    for i in range(NUMCHIPS):
      xline = int((os.popen("echo ${OFFSETX} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      lab.axvline(x=xline, color='k', lw=0.5)
      xline2 = xline + CHIPXMAX
      lab.axvline(x=xline2, color='k', lw=0.5)
      yline = int((os.popen("echo ${OFFSETY} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      plt.axhline(y=yline, color='k', lw=0.5)
      yline2 = yline + CHIPYMAX
      lab.axhline(y=yline2, color='k', lw=0.5)


def calculatingeps(i):
  print("Start calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  A = float(((os.popen("cat " + path + "chip_all.dat | grep A | awk '{print $3}'").readlines())[0]).strip())
  B = float(((os.popen("cat " + path + "chip_all.dat | grep B | awk '{print $3}'").readlines())[0]).strip())
  C = float(((os.popen("cat " + path + "chip_all.dat | grep C | awk '{print $3}'").readlines())[0]).strip())
  D = float(((os.popen("cat " + path + "chip_all.dat | grep D | awk '{print $3}'").readlines())[0]).strip())
  E = float(((os.popen("cat " + path + "chip_all.dat | grep E | awk '{print $3}'").readlines())[0]).strip())
  FCHIP = float(((os.popen("cat " + path + "chip_all.dat | grep -m1 F" + str(i+1) + " | awk '{print $3}'").readlines())[0]).strip())

  x = 2.0*(np.arange(0, int(CHIPXMAX), 1)+offset[i][1])/PIXXMAX
  xred = 2.0*(np.arange(0, int(CHIPXMAX), 10)+offset[i][1])/PIXXMAX
  X = np.arange(0, int(CHIPXMAX), 10)+offset[i][1]

  y = 2.0*(np.arange(0, int(CHIPYMAX), 1)+offset[i][2])/PIXYMAX
  yred = 2.0*(np.arange(0, int(CHIPYMAX), 10)+offset[i][2])/PIXYMAX
  Y = np.arange(0, int(CHIPYMAX), 10)+offset[i][2]

  xx, yy = np.meshgrid(x, y)
  xxred, yyred = np.meshgrid(xred, yred)
  XX, YY = np.meshgrid(X, Y)

  epswZP = A * xx**2 + B * yy**2 + C * xx * yy + D * xx + E * yy + FCHIP
  epswZPred = A * xxred**2 + B * yyred**2 + C * xxred * yyred + D * xxred + E * yyred + FCHIP
  epswoZP = epswZPred - FCHIP
  flux = pow(10,(epswZP/(-2.5)))
  
  CHIP = np.array(len(XX.flatten())*[FCHIP])

  print("Finish calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  print("Start creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  hdu = pyfits.PrimaryHDU(data=flux)
  hdu.writeto(path + 'chip_%i.fits' %(i+1))
  print("Finish creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  return(XX, YY, epswoZP, CHIP)


def plot(arg):
  if (arg == 'resfit'):
    lab.clf()
    fig = plt.figure()

    #Plotting detected objects from the side of the x-axis (fitted)
    ax1 = fig.add_subplot(2,2,1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=True)
    ax1.plot(XXALL,epswoZPALL,'k.')
    appearance('','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax1, 'sub')

    ax3 = fig.add_subplot(2,2,3)
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    ax3.plot(XXALL,FCHIPS,'k.')
    appearance('Xpos','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax3, 'sub')

    #Plotting detected objects from the side of the y-axis (fitted)
    ax2 = fig.add_subplot(2,2,2)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.plot(YYALL,epswoZPALL,'k.')
    appearance('','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax2, 'sub')

    ax4 = fig.add_subplot(2,2,4)
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    ax4.plot(YYALL,FCHIPS,'k.')
    appearance('Ypos','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax4, 'sub')

    lab.savefig(path + 'residuals_fitfunction.png')

  elif (arg == 'cam'):
    lab.clf()
    fig = plt.figure()

    xi = lab.linspace(min(XXALL), max(XXALL))
    yi = lab.linspace(min(YYALL), max(YYALL))

    #Plotting used objects
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(d[:,12],d[:,13],'k,')
    appearance('','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')

    #Plotting applied contour fit
    ax2 = fig.add_subplot(2, 2, 2)
    zi = lab.griddata(XXALL, YYALL, epswZPALL, xi, yi)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.contourf(xi, yi, zi, 10)
    axbar = ax2.contourf(xi, yi, zi, 10)
    fig.colorbar(axbar)
    appearance('','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')

    ax3 = fig.add_subplot(2, 2, 3)
    zi = lab.griddata(XXALL, YYALL, epswoZPALL, xi, yi)
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    ax3.contourf(xi, yi, zi, 10)
    axbar = ax3.contourf(xi, yi, zi, 10)
    fig.colorbar(axbar)
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')

    ax4 = fig.add_subplot(2, 2, 4)
    zi = lab.griddata(XXALL, YYALL, FCHIPS, xi, yi)
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    ax4.contourf(xi, yi, zi, 10)
    axbar = ax4.contourf(xi, yi, zi, 10)
    fig.colorbar(axbar)
    appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')
    
    lab.savefig(path + 'camera.png')
  
  elif (arg == 'mag'):
    lab.clf()
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(d[:,14],d[:,6],'k,')
    appearance('','Residual','', '','', 'grid', 'nocamgrid', ax1, 'sub')

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(d[:,9],d[:,10],'k,')
    appearance('Mag','Residual','', '','', 'grid', 'nocamgrid', ax2, 'sub')
    
    lab.savefig(path + 'mag_dependency.png')

  elif (arg == 'res'):
    lab.clf()
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(d[:,12],d[:,6],'k,')
    appearance('Xpos','Residual','Before fitting', [CAMXMIN, CAMXMAX], [LRL, URL], 'grid', 'nocamgrid', ax1, 'sub')

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(d[:,12],d[:,10],'k,')
    appearance('Xpos','','After fitting', [CAMXMIN, CAMXMAX], [LRL, URL], 'grid', 'nocamgrid', ax2, 'sub')

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(d[:,13],d[:,6],'k,')
    appearance('Ypos','Residual','', [CAMYMIN, CAMYMAX], [LRL, URL], 'grid', 'nocamgrid', ax3, 'sub')

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(d[:,13],d[:,10],'k,')
    appearance('Ypos','','', [CAMYMIN, CAMYMAX], [LRL, URL], 'grid', 'nocamgrid', ax4, 'sub')

    lab.savefig(path + 'residuals.png')
    
  elif (arg == 'histo'):
    lab.clf()
    fig = plt.figure()
    BIN = int((abs(LRL)+abs(URL))/0.01/2.0)+3

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.hist(d[:,10],bins=BIN, range=(LRL,URL))
    appearance('Residual','','After fitting', [LRL, URL], '', 'grid', 'nocamgrid', ax2, 'sub')

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.hist(d[:,6],bins=BIN, range=(LRL,URL))
    appearance('Residual','Number of objects','Before fitting', [LRL, URL], ax2.set_ylim(), 'grid', 'nocamgrid', ax1, 'sub')

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.hist(d[:,10],bins=41, range=(-0.105, 0.105))
    appearance('Residual','','', [-0.1,0.1], '', 'grid', 'nocamgrid', ax4, 'sub')

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.hist(d[:,6],bins=41, range=(-0.105, 0.105))
    appearance('Residual','Number of objects','', [-0.1,0.1], ax4.set_ylim(), 'grid', 'nocamgrid', ax3, 'sub')

    lab.savefig(path + 'histograms.png')
    
  elif (arg == 'contour' ):
    lab.clf()
    fig = plt.figure()
    steps = (abs(LRL)+abs(URL))/10.0
    levels_limit = [np.round(-abs(LRL),2), -0.05, 0.05, np.round(abs(URL),2)]
    xi = lab.linspace(min(d[:,12]), max(d[:,12]))
    yi = lab.linspace(min(d[:,13]), max(d[:,13]))

    ax1 = fig.add_subplot(2, 2, 1)
    zi = lab.griddata(d[:,12], d[:,13], d[:,6], xi, yi)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=True)
    axbar = ax1.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2),color=k)
    fig.colorbar(axbar)
    ax1.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
    appearance('','Ypos','Before fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')

    ax2 = fig.add_subplot(2, 2, 2)
    zi = lab.griddata(d[:,12], d[:,13], d[:,10], xi, yi)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    axbar = ax2.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
    fig.colorbar(axbar)
    ax2.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
    appearance('','','After fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')

    ax3 = fig.add_subplot(2, 2, 3)
    zi = lab.griddata(d[:,12], d[:,13], d[:,6], xi, yi)
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    axbar = ax3.contourf(xi, yi, zi,levels=levels_limit)
    fig.colorbar(axbar)
    ax3.contourf(xi, yi, zi,levels=levels_limit)
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')

    ax4 = fig.add_subplot(2, 2, 4)
    zi = lab.griddata(d[:,12], d[:,13], d[:,10], xi, yi)
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    axbar = ax4.contourf(xi, yi, zi,levels=levels_limit)
    fig.colorbar(axbar)
    ax4.contourf(xi, yi, zi,levels=levels_limit)
    appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')

    lab.savefig(path + 'contourplots.png')
    
    
    
if __name__ == '__main__':
  global offset
  global path
  global CHIPXMAX
  global ROWMAX
  global CHIPYMAX
  global COLUMNMAX
  global PIXXMAX
  global PIXYMAX
  global XXALL
  global YYALL
  global epswoZPALL
  global epswZPALL
  global FCHIPS
  global CAMXMIN
  global CAMXMAX
  global CAMYMIN
  global CAMYMAX
  global NUMCHIPS
  global OFFSETX
  global OFFSETY
  global LRL
  global URL

  
  #Reading command line arguments
  path = sys.argv[1]
  LRL = float(sys.argv[2])
  URL = float(sys.argv[3])

  #Reading chip geometry from config file
  NUMCHIPS = int((os.popen("echo ${NCHIPS} | awk '{print $1}'").readlines())[0])
  PIXXMIN = 0
  PIXXMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $1}'").readlines())[0]) * int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
  PIXYMIN = 0
  PIXYMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $2}'").readlines())[0]) * int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])
  ROWMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $1}'").readlines())[0])
  COLUMNMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $2}'").readlines())[0])
  CHIPXMAX = float((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
  CHIPYMAX = float((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])
  CAMXMIN = int((os.popen("echo ${OFFSETX} | awk '{print $1}'").readlines())[0])
  CAMXMAX = int((os.popen("echo ${OFFSETX} | awk '{print $NF}'").readlines())[0]) + int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
  CAMYMIN = int((os.popen("echo ${OFFSETY} | awk '{print $1}'").readlines())[0])
  CAMYMAX = int((os.popen("echo ${OFFSETY} | awk '{print $NF}'").readlines())[0]) + int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])


  #Reading chip offsets from config file
  offset=([])
  for k in range(NUMCHIPS):
    a = int((os.popen("echo ${OFFSETX} | awk '{print $%i}'" %(k+1)).readlines())[0])
    b = int((os.popen("echo ${OFFSETY} | awk '{print $%i}'" %(k+1)).readlines())[0])
    offset = np.append(offset,[k,a,b])
  offset = offset.reshape((-1,3))


  #get the number of CPUs / cores
  totalcpus = multiprocessing.cpu_count()
  maxcpus = int((os.popen("echo ${NPARA}").readlines())[0])
  usedcpus = maxcpus

  if NUMCHIPS < usedcpus:
    usedcpus = NUMCHIPS

  print("Start calculating with " + str(usedcpus) + "/" + str(totalcpus) + " CPUs...")

  # Initialise process pool:
  pool = multiprocessing.Pool(usedcpus)

  # execute the conversion with a 'pool-map' command:
  catlist = []
  for h in range(NUMCHIPS):
	  catlist.append(h)

  ALL = pool.map(calculatingeps, catlist)


  print("Start plotting...")


  
  XXALL = np.array([])
  YYALL = np.array([])
  epswoZPALL = np.array([])
  FCHIPS = np.array([])
  epswZPALL = np.array([])

  for i in range(NUMCHIPS):
    XXALL = np.append(XXALL,ALL[i][0])
    YYALL = np.append(YYALL,ALL[i][1])
    epswoZPALL = np.append(epswoZPALL,ALL[i][2])
    FCHIPS = np.append(FCHIPS, ALL[i][3])
  epswZPALL = epswoZPALL + FCHIPS

  c = np.array([])
  for k in range(NUMCHIPS):
    c = np.append(c,np.fromfile(path + "chip_%i.csv" %(k+1), sep="\t"))
  d = c.reshape((-1,15))

  pool = multiprocessing.Pool(usedcpus)
  plotlist = ['resfit', 'cam', 'mag', 'res', 'histo', 'contour']
  pool.map(plot, plotlist)