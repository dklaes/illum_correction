#!/usr/bin/python
# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

#Importing packages
from __future__ import division, print_function
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
import os
import sys

#Reading command line arguments
path = sys.argv[1]
LRL = float(sys.argv[2])
URL = float(sys.argv[3])

#Reading chip geometry from config file
NUMCHIPS = int((os.popen("echo ${NCHIPS} | awk '{print $1}'").readlines())[0])
ROWMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $1}'").readlines())[0])
COLUMNMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $2}'").readlines())[0])
CHIPXMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
CHIPYMAX = int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])
CAMXMIN = int((os.popen("echo ${OFFSETX} | awk '{print $1}'").readlines())[0])
CAMXMAX = int((os.popen("echo ${OFFSETX} | awk '{print $NF}'").readlines())[0]) + int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $3}'").readlines())[0])
CAMYMIN = int((os.popen("echo ${OFFSETY} | awk '{print $1}'").readlines())[0])
CAMYMAX = int((os.popen("echo ${OFFSETY} | awk '{print $NF}'").readlines())[0]) + int((os.popen("echo ${CHIPGEOMETRY} | awk '{print $4}'").readlines())[0])

def appearance(xlabel,ylabel,title, xlimits, ylimits, grid, bar, camgrid, plt, sub):
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
    if (bar == 'colorbar'):
	fig.colorbar(axbar)

  if (camgrid == 'camgrid'):
    global OFFSETX
    global OFFSETY
    global CHIPXMAX
    global CHIPYMAX
    for i in range(NUMCHIPS):
      xline = int((os.popen("echo ${OFFSETX} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      lab.axvline(x=xline, color='k', lw=0.5)
      xline2 = xline + CHIPXMAX
      lab.axvline(x=xline2, color='k', lw=0.5)
      yline = int((os.popen("echo ${OFFSETY} | awk '{print $" + str(i+1) + "}'").readlines())[0])
      plt.axhline(y=yline, color='k', lw=0.5)
      yline2 = yline + CHIPYMAX
      lab.axhline(y=yline2, color='k', lw=0.5)

print("Plotting fitted residuals from X- and Y-direction and residuals vs. airmass...")

#Plotting detected objects from the side of the x-axis
a = np.array([])
for k in range(NUMCHIPS):
  a = np.append(a,np.fromfile(path + "chip_%i.csv" %(k+1), sep="\t"))
b = a.reshape((-1,14))

####Plotting airmass (x-axis) vs. residuals (y-axis) after fitting
###lab.clf()
###lab.plot(b[:,11],b[:,10],'k,')
###appearance('Airmass','Residual','Residuals', '', '', 'grid', 'nocolorbar', 'nocamgrid', plt, 'nosub')
###lab.savefig(path + 'Airmass_afterFitting.png')

###lab.clf()
###lab.plot(b[:,11],b[:,6],'k,')
###appearance('Airmass','Residual','Residuals', '', '', 'grid', 'nocolorbar', 'nocamgrid', plt, 'nosub')
###lab.savefig(path + 'Airmass_beforeFitting.png')



lab.clf()
fig = plt.figure()

ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(b[:,12],b[:,6],'k,')
appearance('Xpos','Residual','Before fitting', [CAMXMIN, CAMXMAX], [LRL, URL], 'grid', 'nocolorbar', 'nocamgrid', ax1, 'sub')

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(b[:,12],b[:,10],'k,')
appearance('Xpos','','After fitting', [CAMXMIN, CAMXMAX], [LRL, URL], 'grid', 'nocolorbar', 'nocamgrid', ax2, 'sub')

ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(b[:,13],b[:,6],'k,')
appearance('Ypos','Residual','', [CAMYMIN, CAMYMAX], [LRL, URL], 'grid', 'nocolorbar', 'nocamgrid', ax3, 'sub')

ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(b[:,13],b[:,10],'k,')
appearance('Ypos','','', [CAMYMIN, CAMYMAX], [LRL, URL], 'grid', 'nocolorbar', 'nocamgrid', ax4, 'sub')

lab.savefig(path + 'residuals.png')



lab.clf()
fig = plt.figure()
BIN = int((abs(LRL)+abs(URL))/0.01/2.0)+3

ax2 = fig.add_subplot(2, 2, 2)
ax2.hist(b[:,10],bins=BIN, range=(LRL,URL))
appearance('Residual','','After fitting', [LRL, URL], '', 'grid', 'nocolorbar', 'nocamgrid', ax2, 'sub')

ax1 = fig.add_subplot(2, 2, 1)
ax1.hist(b[:,6],bins=BIN, range=(LRL,URL))
appearance('Residual','Number of objects','Before fitting', [LRL, URL], ax2.set_ylim(), 'grid', 'nocolorbar', 'nocamgrid', ax1, 'sub')

ax4 = fig.add_subplot(2, 2, 4)
ax4.hist(b[:,10],bins=41, range=(-0.105, 0.105))
appearance('Residual','','', [-0.1,0.1], '', 'grid', 'nocolorbar', 'nocamgrid', ax4, 'sub')

ax3 = fig.add_subplot(2, 2, 3)
ax3.hist(b[:,6],bins=41, range=(-0.105, 0.105))
appearance('Residual','Number of objects','', [-0.1,0.1], ax4.set_ylim(), 'grid', 'nocolorbar', 'nocamgrid', ax3, 'sub')

lab.savefig(path + 'histograms.png')



lab.clf()
fig = plt.figure()
steps = (abs(LRL)+abs(URL))/10.0
levels_limit = [np.round(-abs(LRL),2), -0.05, 0.05, np.round(abs(URL),2)]

ax1 = fig.add_subplot(2, 2, 1)
xi = lab.linspace(min(b[:,12]), max(b[:,12]))
yi = lab.linspace(min(b[:,13]), max(b[:,13]))
zi = lab.griddata(b[:,12], b[:,13], b[:,6], xi, yi)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=True)
ax1.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
axbar = ax1.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2),color=k)
appearance('','Ypos','Before fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax1, 'sub')

ax2 = fig.add_subplot(2, 2, 2)
xi = lab.linspace(min(b[:,12]), max(b[:,12]))
yi = lab.linspace(min(b[:,13]), max(b[:,13]))
zi = lab.griddata(b[:,12], b[:,13], b[:,10], xi, yi)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
axbar = ax2.contourf(xi, yi, zi,levels=np.round(np.arange(LRL, URL+(steps/2.0), steps),2))
appearance('','','After fitting', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax2, 'sub')

ax3 = fig.add_subplot(2, 2, 3)
xi = lab.linspace(min(b[:,12]), max(b[:,12]))
yi = lab.linspace(min(b[:,13]), max(b[:,13]))
zi = lab.griddata(b[:,12], b[:,13], b[:,6], xi, yi)
ax3.contourf(xi, yi, zi,levels=levels_limit)
plt.setp(ax3.get_xticklabels(), visible=True)
plt.setp(ax3.get_yticklabels(), visible=True)
axbar = ax3.contourf(xi, yi, zi,levels=levels_limit)
appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax3, 'sub')

ax4 = fig.add_subplot(2, 2, 4)
xi = lab.linspace(min(b[:,12]), max(b[:,12]))
yi = lab.linspace(min(b[:,13]), max(b[:,13]))
zi = lab.griddata(b[:,12], b[:,13], b[:,10], xi, yi)
ax4.contourf(xi, yi, zi,levels=levels_limit)
plt.setp(ax4.get_xticklabels(), visible=True)
plt.setp(ax4.get_yticklabels(), visible=False)
axbar = ax4.contourf(xi, yi, zi,levels=levels_limit)
appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax4, 'sub')

lab.savefig(path + 'contourplots.png')