# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

# ----------------------------------------------------------------
# File Name:           illum_correction_contourplot_fitfunction.py
# Author:              Dominik Klaes (dklaes@astro.uni-bonn.de)
# Last modified on:    31.05.2013
# Version:             V1.2
# Description:         Create FITS-correction files and corresponding plots
# ----------------------------------------------------------------

# Changes from V1.0 to V1.1
# - more comments
# - Checkplot magnitude / residual dependency added
# - Now running also without X-server
# - old "illum_correction_plot_fitted.py" included and restructed to multiprocessing

# Changes from V1.1 to V1.2
# - scatter plot added

# Changes from V1.2 to V1.3
# - Including single plots
# - Changed condition for colorbars of contour plots from sigma < 0.05 to sigma < 0.06
# - Closed every plot directly after savefig to save memory

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
import getopt
import multiprocessing
import pyfits
import ldac
from matplotlib import cm


def appearance(xlabel, ylabel, title, xlimits, ylimits, grid, camgrid, plt, sub):
  if (grid == 'grid'):
    plt.grid(True)		# Adding grid
  elif (grid == 'nogrid'):
    plt.grid(False)		# Remove grid
  if (sub == 'nosub'):		# For one plot (non-subplot)
    plt.xlabel(xlabel)		# x label
    plt.ylabel(ylabel)		# y label
    plt.title(title)		# Adding a title
    if (xlimits != ''):		# Giving xrange, if empty then auto
	plt.xlim(xlimits)
    if (ylimits != ''):		# Giving yrange, if empty then auto
	plt.ylim(ylimits)
    if (bar == 'colorbar'):
	plt.colorbar()		# Adding colorbar
  elif (sub == 'sub'):		# For subplots
    plt.set_xlabel(xlabel)	# x label
    plt.set_ylabel(ylabel)	# y label
    plt.set_title(title)	# Adding a title
    if (xlimits != ''):		# Giving xrange, if empty then auto
	plt.set_xlim(xlimits)
    if (ylimits != ''):		# Giving yrange, if empty then auto
	plt.set_ylim(ylimits)

  if (camgrid == 'camgrid'):	# Adding a camera grid, one line for each border of each chip
    for i in range(NUMCHIPS):
      xline = int((os.popen("echo ${OFFSETX} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      lab.axvline(x=xline, color='k', lw=0.5)
      xline2 = xline + CHIPXMAX
      lab.axvline(x=xline2, color='k', lw=0.5)
      yline = int((os.popen("echo ${OFFSETY} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      plt.axhline(y=yline, color='k', lw=0.5)
      yline2 = yline + CHIPYMAX
      lab.axhline(y=yline2, color='k', lw=0.5)


def calculatingeps(input):
  i = input[0]
  # Getting the prefactors.
  A = input[1]
  B = input[2]
  C = input[3]
  D = input[4]
  E = input[5]
  FCHIP = input[6]
  print("Start calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")  

  # Creating arrays for x (resized chip coordinates (for numerical reasons)), xred (reduced and resized chip
  # coordinates (for plotting and numerical reasons)) and X (reduced (for plotting)). Having xred and X seperated doesn't
  # cause the problem of transforming the prefactors to the other coordinate system!
  x = 2.0*(np.arange(0, int(CHIPXMAX), 1)+offset[i][1])/PIXXMAX
  xred = 2.0*(np.arange(0, int(CHIPXMAX), 10)+offset[i][1])/PIXXMAX
  X = np.arange(0, int(CHIPXMAX), 10)+offset[i][1]

  # Now the same for y, yred and Y.
  y = 2.0*(np.arange(0, int(CHIPYMAX), 1)+offset[i][2])/PIXYMAX
  yred = 2.0*(np.arange(0, int(CHIPYMAX), 10)+offset[i][2])/PIXYMAX
  Y = np.arange(0, int(CHIPYMAX), 10)+offset[i][2]

  # Creating the corresponding meshgrids.
  xx, yy = np.meshgrid(x, y)
  xxred, yyred = np.meshgrid(xred, yred)
  XX, YY = np.meshgrid(X, Y)

  # Calculating the position dependend residuals.
  # - epswZP:		residuals with chip zeropoints in resized chip coordinates,
  #			for FITS correction file
  # - flux:		multiplicative factor to correct illumination effects
  #			for FITS correction file
  # - epswZPred:	residuals with chip zeropoints in resized chip coordinates (reduced data set),
  #			for plotting
  # - epswoZP:		residuals without chip zeropoints in resized chip coordinates (reduced data set),
  #			for plotting
  epswZP = A * xx**2 + B * yy**2 + C * xx * yy + D * xx + E * yy + FCHIP
  flux = pow(10,(epswZP/(-2.5)))
  epswZPred = A * xxred**2 + B * yyred**2 + C * xxred * yyred + D * xxred + E * yyred + FCHIP
  epswoZP = epswZPred - FCHIP
  
  # Creating an array containing the chip zeropoint offset for each data point (for plotting reasons).
  CHIP = np.array(len(XX.flatten())*[FCHIP])

  print("Finish calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  
  # Writing the correction FITS files.
  print("Start creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  hdu = pyfits.PrimaryHDU(data=flux)
  hdu.writeto(path + 'chip_%i.fits' %(i+1))
  print("Finish creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  
  # Returning data for plotting.
  return(XX, YY, epswoZP, CHIP)


def plot(arg):
  if (arg == 'resfit'):
    # Plotting how the fitting function actually looks like. Plots are:
    # 1:	Upper left:	Fitting function without zeropoints, from x-axis
    # 2:	Upper right:	Fitting function without zeropoints, from y-axis
    # 3:	Lower left:	zeropoints, from x-axis
    # 4:	Lower right:	zeropoints, from y-axis
    
    lab.clf()
    fig = plt.figure()

    # 1
    ax1 = fig.add_subplot(2,2,1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=True)
    ax1.plot(XXALL,epswoZPALL,'k.')
    appearance('','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax1, 'sub')

    # 2
    ax2 = fig.add_subplot(2,2,2)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.plot(YYALL,epswoZPALL,'k.')
    appearance('','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax2, 'sub')
    
    # 3
    ax3 = fig.add_subplot(2,2,3)
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    ax3.plot(XXALL,FCHIPS,'k.')
    appearance('Xpos','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax3, 'sub')

    # 4
    ax4 = fig.add_subplot(2,2,4)
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    ax4.plot(YYALL,FCHIPS,'k.')
    appearance('Ypos','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax4, 'sub')

    lab.savefig(path + 'residuals_fitfunction.png')
    plt.close()


    # Plot now all plots individually in single files.

    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(XXALL,epswoZPALL,'k.')
    appearance('Xpos','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_fitfunction_1.png')
    plt.close()
    
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1,1,1)
    ax2.plot(YYALL,epswoZPALL,'k.')
    appearance('Ypos','Residual','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_fitfunction_2.png')
    plt.close()
    
    # 3
    lab.clf()
    fig = plt.figure()
    ax3 = fig.add_subplot(1,1,1)
    ax3.plot(XXALL,FCHIPS,'k.')
    appearance('Xpos','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocamgrid', ax3, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_fitfunction_3.png')
    plt.close()
    
    # 4
    lab.clf()
    fig = plt.figure()
    ax4 = fig.add_subplot(1,1,1)
    ax4.plot(YYALL,FCHIPS,'k.')
    appearance('Ypos','Residual','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocamgrid', ax4, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_fitfunction_4.png')
    plt.close()


  elif (arg == 'cam'):
    # Plotting camera related stuff. Plots are:
    # 1:	Upper left:	Used objects are being plotted.
    # 2:	Upper right:	Correction as contour plot (with zeropoints)
    # 3:	Lower left:	Correction as contour plot (without zeropoints)
    # 4:	Lower right:	Correction as contour plot (only zeropoints)
    
    lab.clf()
    fig = plt.figure()

    xi = lab.linspace(min(XXALL), max(XXALL))
    yi = lab.linspace(min(YYALL), max(YYALL))

    # 1
    ax1 = fig.add_subplot(2, 2, 1, aspect='equal')
    ax1.plot(data['Xpos_global'],data['Ypos_global'],'k,')
    appearance('','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')

    # 2
    ax2 = fig.add_subplot(2, 2, 2, aspect='equal')
    zi2 = lab.griddata(XXALL, YYALL, epswZPALL, xi, yi)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    CS2 = ax2.contourf(xi, yi, zi2, 10)
    cbar2 = plt.colorbar(CS2, format='%0.2f')
    cbar2.ax.set_ylabel('Residual in mag')
    appearance('','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')

    # 3
    ax3 = fig.add_subplot(2, 2, 3, aspect='equal')
    zi3 = lab.griddata(XXALL, YYALL, epswoZPALL, xi, yi)
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    CS3 = ax3.contourf(xi, yi, zi3, 10)
    cbar3 = plt.colorbar(CS3, format='%0.2f')
    cbar3.ax.set_ylabel('Residual in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')

    # 4
    ax4 = fig.add_subplot(2, 2, 4, aspect='equal')
    zi4 = lab.griddata(XXALL, YYALL, FCHIPS, xi, yi)
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    CS4 = ax4.contourf(xi, yi, zi4, 10)
    cbar4 = plt.colorbar(CS4, format='%0.2f')
    cbar4.ax.set_ylabel('Residual in mag')
    appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')

    lab.savefig(path + 'camera.png')
    plt.close()
    
    
    # Plot now all plots individually in single files.
    
    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    ax1.plot(data['Xpos_global'],data['Ypos_global'],'k,')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'camera_1.png')
    plt.close()
    
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1, 1, 1, aspect='equal')
    CS2 = ax2.contourf(xi, yi, zi2, 10)
    cbar2 = plt.colorbar(CS2, format='%0.2f')
    cbar2.ax.set_ylabel('Residual in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'camera_2.png')
    plt.close()
    
    # 3
    lab.clf()
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1, aspect='equal')
    CS3 = ax3.contourf(xi, yi, zi3, 10)
    cbar3 = plt.colorbar(CS3, format='%0.2f')
    cbar3.ax.set_ylabel('Residual in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')
    lab.savefig(path + 'single_plots/' + 'camera_3.png')
    plt.close()
    
    # 4
    lab.clf()
    fig = plt.figure()
    ax4 = fig.add_subplot(1, 1, 1, aspect='equal')
    CS4 = ax4.contourf(xi, yi, zi4, 10)
    cbar4 = plt.colorbar(CS4, format='%0.2f')
    cbar4.ax.set_ylabel('Residual in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')
    lab.savefig(path + 'single_plots/' + 'camera_4.png')
    plt.close()
    
  
  elif (arg == 'mag'):
    # Plotting mag-residuals dependency before and after fitting. Plots are:
    # 1:	Upper:	Before fitting
    # 2:	Lower:	After fitting
    
    lab.clf()
    fig = plt.figure()

    # 1
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(data['MagZP'],data['Residual'],'k,')
    appearance('','Residual','', '','', 'grid', 'nocamgrid', ax1, 'sub')

    # 2
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(data['Mag_fitted'],data['Residual_fitted'],'k,')
    appearance('Mag','Residual','', '','', 'grid', 'nocamgrid', ax2, 'sub')
    
    lab.savefig(path + 'mag_dependency.png')
    plt.close()
    
    
    # Plot now all plots individually in single files.
    
    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(data['MagZP'],data['Residual'],'k,')
    appearance('Mag','Residual','', '','', 'grid', 'nocamgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'mag_dependency_1.png')
    plt.close()
    
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1, 1, 1)
    ax2.plot(data['Mag_fitted'],data['Residual_fitted'],'k,')
    appearance('Mag','Residual','', '','', 'grid', 'nocamgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'mag_dependency_2.png')
    plt.close()
    

  elif (arg == 'res'):
    # Residuals shown before and after fitting, from x- and y-axis. Plots are:
    # 1:	Upper left:	Before fitting, from x-axis
    # 2:	Upper right:	After fitting, from x-axis
    # 3:	Lower left:	Before fitting, from y-axis
    # 4:	Lower right:	After fitting, from y-axis

    #datamax = np.amax((np.fabs(d[:,6]), np.fabs(d[:,10])))
    datamax = np.amax((np.fabs(data['Residual']), np.fabs(data['Residual_fitted'])))
    
    lab.clf()
    fig = plt.figure()

    # 1
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(data['Xpos_global'],data['Residual'],'k,')
    appearance('Xpos','Residual','Before fitting', [CAMXMIN, CAMXMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax1, 'sub')

    # 2
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(data['Xpos_global'],data['Residual_fitted'],'k,')
    appearance('Xpos','','After fitting', [CAMXMIN, CAMXMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax2, 'sub')

    # 3
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(data['Ypos_global'],data['Residual'],'k,')
    appearance('Ypos','Residual','', [CAMYMIN, CAMYMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax3, 'sub')

    # 4
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(data['Ypos_global'],data['Residual_fitted'],'k,')
    appearance('Ypos','','', [CAMYMIN, CAMYMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax4, 'sub')

    lab.savefig(path + 'residuals.png')
    plt.close()
    
    
    # Plot now all plots individually in single files.
    
    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(data['Xpos_global'],data['Residual'],'k,')
    appearance('Xpos','Residual','Before fitting', [CAMXMIN, CAMXMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_1.png')
    plt.close()
        
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1, 1, 1)
    ax2.plot(data['Xpos_global'],data['Residual_fitted'],'k,')
    appearance('Xpos','Residual','After fitting', [CAMXMIN, CAMXMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_2.png')
    plt.close()
    
    # 3
    lab.clf()
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1)
    ax3.plot(data['Ypos_global'],data['Residual'],'k,')
    appearance('Ypos','Residual','Before fitting', [CAMYMIN, CAMYMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax3, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_3.png')
    plt.close()
    
    # 4
    lab.clf()
    fig = plt.figure()
    ax4 = fig.add_subplot(1, 1, 1)
    ax4.plot(data['Ypos_global'],data['Residual_fitted'],'k,')
    appearance('Ypos','Residual','After fitting', [CAMYMIN, CAMYMAX], [-datamax, datamax], 'grid', 'nocamgrid', ax4, 'sub')
    lab.savefig(path + 'single_plots/' + 'residuals_4.png')
    plt.close()
    
    
  elif (arg == 'histo'):
    # Histograms of number of objects in magnitude bins, total and zoomed in. Plots are:
    # 1:	Upper left:	Before fitting, not zoomed
    # 2:	Upper right:	After fitting, not zoomed
    # 3:	Lower left:	Before fitting, zoomed into -0.1 to 0.1
    # 4:	Lower right:	After fitting, toomed into -0.1 to 0.1
    
    lab.clf()
    fig = plt.figure()
    BIN = int((abs(-sigma)+abs(sigma))/0.01/2.0)+3

    # Optimized so that y ranges (number of objects) are the same in both plots, the second one has, due to fitting
    # a higher number so we have to take this (automatically) set (y) limits.
    # 2
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.hist(data['Residual_fitted'],bins=BIN, range=(-sigma,sigma))
    appearance('Residual','','After fitting', [-sigma, sigma], '', 'grid', 'nocamgrid', ax2, 'sub')
    
    # 1
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.hist(data['Residual'],bins=BIN, range=(-sigma,sigma))
    appearance('Residual','Number of objects','Before fitting', [-sigma, sigma], ax2.set_ylim(), 'grid', 'nocamgrid', ax1, 'sub')
    
    # 4
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.hist(data['Residual_fitted'],bins=41, range=(-0.105, 0.105))
    appearance('Residual','','', [-0.1,0.1], '', 'grid', 'nocamgrid', ax4, 'sub')

    # 3
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.hist(data['Residual'],bins=41, range=(-0.105, 0.105))
    appearance('Residual','Number of objects','', [-0.1,0.1], ax4.set_ylim(), 'grid', 'nocamgrid', ax3, 'sub')
    
    lab.savefig(path + 'histograms.png')
    plt.close()
    
    
    # Plot now all plots individually in single files.
    
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1, 1, 1)
    ax2.hist(data['Residual_fitted'],bins=BIN, range=(-sigma,sigma))
    appearance('Residual','Number of objects','After fitting', [-sigma, sigma], '', 'grid', 'nocamgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'histograms_2.png')
    plt.close()
    
    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.hist(data['Residual'],bins=BIN, range=(-sigma,sigma))
    appearance('Residual','Number of objects','Before fitting', [-sigma, sigma], ax2.set_ylim(), 'grid', 'nocamgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'histograms_1.png')
    plt.close()
    
    # 4
    lab.clf()
    fig = plt.figure()
    ax4 = fig.add_subplot(1, 1, 1)
    ax4.hist(data['Residual_fitted'],bins=41, range=(-0.105, 0.105))
    appearance('Residual','Number of objects','After fitting', [-0.1,0.1], '', 'grid', 'nocamgrid', ax4, 'sub')
    lab.savefig(path + 'single_plots/' + 'histograms_4.png')
    plt.close()
    
    # 3
    lab.clf()
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1)
    ax3.hist(data['Residual'],bins=41, range=(-0.105, 0.105))
    appearance('Residual','Number of objects','Before fitting', [-0.1,0.1], ax4.set_ylim(), 'grid', 'nocamgrid', ax3, 'sub')
    lab.savefig(path + 'single_plots/' + 'histograms_3.png')
    plt.close()
    
    
  elif (arg == 'contour'):
    # Contourplots of residuals before and after plotting with two different colourbar ranges. Plots are:
    # 1:	Upper left:	Before fitting, auto ranges
    # 2:	Upper right:	After fitting, auto ranges
    # 3:	Lower left:	Before fitting, manuel set ranges to "max < 0.05 < -0.05 < min"
    # 4:	Lower right:	After fitting, manuel set ranges to "max < 0.05 < -0.05 < min"
    
    lab.clf()
    fig = plt.figure()
    steps = (abs(-sigma)+abs(sigma))/10.0
    if (sigma < 0.06):
      levels_limit = [np.round(-abs(minimum),2), -0.05, 0.05, np.round(abs(maximum),2)]
    else:
      levels_limit = [np.round(-abs(-sigma),2), -0.05, 0.05, np.round(abs(sigma),2)]
    xi = lab.linspace(min(data['Xpos_global']), max(data['Ypos_global']))
    yi = lab.linspace(min(data['Xpos_global']), max(data['Ypos_global']))

    # Optimised order so that zi is computed only once.
    # 1
    ax1 = fig.add_subplot(2, 2, 1, aspect='equal')
    zi13 = lab.griddata(data['Xpos_global'], data['Ypos_global'], data['Residual'], xi, yi)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=True)
    axbar = ax1.contourf(xi, yi, zi13,levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2),color=k)
    fig.colorbar(axbar)
    ax1.contourf(xi, yi, zi13, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    appearance('','Ypos','Before fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')

    # 3
    ax3 = fig.add_subplot(2, 2, 3, aspect='equal')
    plt.setp(ax3.get_xticklabels(), visible=True)
    plt.setp(ax3.get_yticklabels(), visible=True)
    axbar = ax3.contourf(xi, yi, zi13,levels=levels_limit)
    fig.colorbar(axbar)
    ax3.contourf(xi, yi, zi13, levels=levels_limit)
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')
    
    
    # Optimised order so that zi is computed only once.
    # 2
    ax2 = fig.add_subplot(2, 2, 2, aspect='equal')
    zi24 = lab.griddata(data['Xpos_global'], data['Ypos_global'], data['Residual_fitted'], xi, yi)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    axbar = ax2.contourf(xi, yi, zi24, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    fig.colorbar(axbar)
    ax2.contourf(xi, yi, zi24, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    appearance('','','After fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')

    # 4
    ax4 = fig.add_subplot(2, 2, 4, aspect='equal')
    plt.setp(ax4.get_xticklabels(), visible=True)
    plt.setp(ax4.get_yticklabels(), visible=False)
    axbar = ax4.contourf(xi, yi, zi24, levels=levels_limit)
    fig.colorbar(axbar)
    ax4.contourf(xi, yi, zi24, levels=levels_limit)
    appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')

    lab.savefig(path + 'contourplots.png')
    plt.close()
    
    
    # Plot now all plots individually in single files.
    
    # 1
    lab.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    axbar = ax1.contourf(xi, yi, zi13, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2),color=k)
    fig.colorbar(axbar)
    ax1.contourf(xi, yi, zi13, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    appearance('Xpos','Ypos','Before fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')
    lab.savefig(path + 'single_plots/' + 'contourplots_1.png')
    plt.close()
    
    # 3
    lab.clf()
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1, aspect='equal')
    axbar = ax3.contourf(xi, yi, zi13,levels=levels_limit)
    fig.colorbar(axbar)
    ax3.contourf(xi, yi, zi13, levels=levels_limit)
    appearance('Xpos','Ypos','Before fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax3, 'sub')
    lab.savefig(path + 'single_plots/' + 'contourplots_3.png')
    plt.close()
    
    # 2
    lab.clf()
    fig = plt.figure()
    ax2 = fig.add_subplot(1, 1, 1, aspect='equal')
    axbar = ax2.contourf(xi, yi, zi24, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    fig.colorbar(axbar)
    ax2.contourf(xi, yi, zi24, levels=np.round(np.arange(-sigma, sigma+(steps/2.0), steps),2))
    appearance('Xpos','Ypos','After fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax2, 'sub')
    lab.savefig(path + 'single_plots/' + 'contourplots_2.png')
    plt.close()
    
    # 4
    lab.clf()
    fig = plt.figure()
    ax4 = fig.add_subplot(1, 1, 1, aspect='equal')
    axbar = ax4.contourf(xi, yi, zi24, levels=levels_limit)
    fig.colorbar(axbar)
    ax4.contourf(xi, yi, zi24, levels=levels_limit)
    appearance('Xpos','Ypos','After fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax4, 'sub')
    lab.savefig(path + 'single_plots/' + 'contourplots_4.png')
    plt.close()
    

  elif (arg == 'scatter'):
    # Scatter plot with location of used objects with color coded depending on residual value before and after fitting.

    lab.clf()
    fig = plt.figure()

    datamax = np.amax((np.fabs(data['Residual']), np.fabs(data['Residual_fitted'])))
    intervall = np.linspace(-datamax, datamax, num=11)

    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    cax = ax1.scatter(data['Xpos_global'], data['Ypos_global'], c=data['Residual'], marker='o', cmap=cm.spectral, vmin=-datamax, vmax=datamax)
    cbar = plt.colorbar(cax, ticks=intervall, format='%0.2f')

    cbar.ax.set_ylabel('Residuum in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')
    
    lab.savefig(path + 'scatter_before_fitting.png')
    lab.savefig(path + 'single_plots/' + 'scatter_before_fitting.png')
    plt.close()


    lab.clf()
    fig = plt.figure()

    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    cax = ax1.scatter(data['Xpos_global'], data['Ypos_global'], c=data['Residual_fitted'], marker='o', cmap=cm.spectral, vmin=-datamax, vmax=datamax)
    cbar = plt.colorbar(cax, ticks=intervall, format='%0.2f')

    cbar.ax.set_ylabel('Residuum in mag')
    appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'camgrid', ax1, 'sub')

    lab.savefig(path + 'scatter_after_fitting.png')
    lab.savefig(path + 'single_plots/' + 'scatter_after_fitting.png')
    plt.close()
    
    
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
  global sigma
  global minimum
  global maxmimum
  global data
  
  
  #Reading command line arguments
  opts, args = getopt.getopt(sys.argv[1:], "i:p:t:e:", ["input=", "path=", "table=", "external="])

  infile = path = table =  external = None
  for o, a in opts:
      if o in ("-i"):
	  infile = a.split()
      elif o in ("-p"):
	  path = a
      elif o in ("-t"):
	  table = a
      elif o in ("-e"):
	  external = a.split()

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

  # Use only as many cores as chips are avaiable
  if NUMCHIPS < usedcpus:
    usedcpus = NUMCHIPS

  print("Start calculating with " + str(usedcpus) + "/" + str(totalcpus) + " CPUs...")

  # Initialise process pool:
  pool = multiprocessing.Pool(usedcpus)

  # Getting the prefactors.
  coefffile = external[0]
  
  coeffs = {}
  f = open(coefffile, 'r')
  for line in f:
       key, val = line.split()[0], float(line.split()[2])
       coeffs[key] = val
  f.close()
  
  A = float(coeffs['A'])
  B = float(coeffs['B'])
  C = float(coeffs['C'])
  D = float(coeffs['D'])
  E = float(coeffs['E'])

  # execute the calculating with a 'pool-map' command:
  catlist = []
  for h in range(NUMCHIPS):
	  FCHIP = float(coeffs['F' + str(h+1)])
	  catlist.append((h, A, B, C, D, E, FCHIP))
  
  ALL = pool.map(calculatingeps, catlist)


  print("Start plotting...")

  # Split the returning "ALL" array from "calculatingeps" into components:
  # - XXALL:		X meshgrid of reduced chip coordinates
  # - YYALL:		Y meshgrid of reduced chip coordinates
  # - epswoZPALL:	residuals without chip zeropoint offsets for reduced chip coordinates
  # - FCHIPS:		Array that contains the zeropoint offsets for reach reduced chip coordinates
  # - wpswZPALL:	residuals with chip zeropoint offsets for reduced chip coordinates
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

  # Importing catalogues.
  data = ldac.LDACCat(infile[0])[table]
  sigma = np.std(data['Residual_fitted'])
  minimum = np.amin(data['Residual_fitted'])
  maximum = np.amax(data['Residual_fitted'])

  # Running plots on multiple cores to save time.
  pool = multiprocessing.Pool(usedcpus)
  plotlist = ['resfit', 'cam', 'mag', 'res', 'histo', 'contour', 'scatter']
  pool.map(plot, plotlist)
