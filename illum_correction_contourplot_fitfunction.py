#!/usr/bin/python
# -*- coding: utf-8 -*-
#scipy-0.11 and numpy-1.6.2 required!

#Importing packages
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
import os
import sys
import multiprocessing
import pyfits

#Reading command line arguments
path = sys.argv[1]

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

def appearance(xlabel, ylabel, title, xlimits, ylimits, grid, bar, camgrid, plt, sub):
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


#Reading chip offsets from config file
offset=([])
for k in range(NUMCHIPS):
  x = int((os.popen("echo ${OFFSETX} | awk '{print $%i}'" %(k+1)).readlines())[0])
  y = int((os.popen("echo ${OFFSETY} | awk '{print $%i}'" %(k+1)).readlines())[0])
  offset = np.append(offset,[k,x,y])
offset = offset.reshape((-1,3))

def calculatingeps(i):
  print("Start calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  global offset
  global path
  global CHIPXMAX
  global ROWMAX
  global CHIPYMAX
  global COLUMNMAX
  global PIXXMAX
  global PIXYMAX
  A = float(((os.popen("cat " + path + "chip_all.dat | grep A | awk '{print $3}'").readlines())[0]).strip())
  B = float(((os.popen("cat " + path + "chip_all.dat | grep B | awk '{print $3}'").readlines())[0]).strip())
  C = float(((os.popen("cat " + path + "chip_all.dat | grep C | awk '{print $3}'").readlines())[0]).strip())
  D = float(((os.popen("cat " + path + "chip_all.dat | grep D | awk '{print $3}'").readlines())[0]).strip())
  E = float(((os.popen("cat " + path + "chip_all.dat | grep E | awk '{print $3}'").readlines())[0]).strip())
  FCHIP = float(((os.popen("cat " + path + "chip_all.dat | grep -m1 F" + str(i+1) + " | awk '{print $3}'").readlines())[0]).strip())
  data = np.zeros((int(CHIPXMAX*CHIPYMAX)))
  e = np.zeros((int(CHIPXMAX*CHIPYMAX*5)))	#*5 for x,y,eps with and without ZP offset and ZP offset itself for each coordinate
  for c in range(0,int(CHIPYMAX),1):
   Y = (2.0*(c+offset[i][2]))/PIXYMAX
   for a in range(0,int(CHIPXMAX),1):
     X = (2.0*(a+offset[i][1]))/PIXXMAX
     epswZP = A * X ** 2 + B * Y ** 2 + C * X * Y + D * X + E * Y + FCHIP
     epswoZP = epswZP - FCHIP
     flux = pow(10,(epswZP/(-2.5)))
     e.itemset((int(c*CHIPXMAX*5+a*5)),a+offset[i][1])
     e.itemset((int(c*CHIPXMAX*5+a*5+1)),c+offset[i][2])
     e.itemset((int(c*CHIPXMAX*5+a*5+2)),epswZP)
     e.itemset((int(c*CHIPXMAX*5+a*5+3)),epswoZP)
     e.itemset((int(c*CHIPXMAX*5+a*5+4)),FCHIP)
     data.itemset((int(c*CHIPXMAX+a)),flux)
  print("Finish calculating data for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  print("Start creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  data = data.reshape((CHIPYMAX,CHIPXMAX))
  hdu = pyfits.PrimaryHDU(data=data)
  hdu.writeto(path + 'chip_%i.fits' %(i+1))
  print("Finish creating FITS-file for chip " + str(i+1) + "/" + str(NUMCHIPS) + "...")
  e = e.flatten()
  e = e.reshape((-1,5))
  return(e[::100])

# get the number of CPUs / cores
totalcpus = multiprocessing.cpu_count()
maxcpus = int((os.popen("echo ${NPARA}").readlines())[0])
usedcpus = maxcpus

if NUMCHIPS < usedcpus:
  usedcpus = NUMCHIPS

os.popen("echo Begin Berechnung: `date` >> ~/dauer.txt")

print("Start calculating with " + str(usedcpus) + "/" + str(totalcpus) + " CPUs...")

# Initialise process pool:
pool = multiprocessing.Pool(usedcpus)

# execute the conversion with a 'pool-map' command:
catlist = []
for h in range(NUMCHIPS):
        catlist.append(h)

b = np.array([])
b = np.append(b, pool.map(calculatingeps, catlist))

os.popen("echo Ende Berechnung: `date` >> ~/dauer.txt")
os.popen("echo Beginn plotten: `date` >> ~/dauer.txt")

print("Start plotting...")
b = b.flatten()
b = b.reshape((-1,5))
indices = np.lexsort(keys = (b[:][1], b[:][0]))
b.take(indices, axis=-1)


#Plotting detected objects from the side of the x-axis (fitted)
lab.clf()
fig = plt.figure()

ax1 = fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=True)
ax1.plot(b[:,0],b[:,3],'k.')
appearance('','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocolorbar', 'nocamgrid', ax1, 'sub')

ax3 = fig.add_subplot(2,2,3)
plt.setp(ax3.get_xticklabels(), visible=True)
plt.setp(ax3.get_yticklabels(), visible=True)
ax3.plot(b[:,0],b[:,4],'k.')
appearance('Xpos','Residual','', [CAMXMIN, CAMXMAX], '', 'grid', 'nocolorbar', 'nocamgrid', ax3, 'sub')


#Plotting detected objects from the side of the y-axis (fitted)
ax2 = fig.add_subplot(2,2,2)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.plot(b[:,1],b[:,3],'k.')
appearance('','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocolorbar', 'nocamgrid', ax2, 'sub')


ax4 = fig.add_subplot(2,2,4)
plt.setp(ax4.get_xticklabels(), visible=True)
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.plot(b[:,1],b[:,4],'k.')
appearance('Ypos','','', [CAMYMIN, CAMYMAX], '', 'grid', 'nocolorbar', 'nocamgrid', ax4, 'sub')
lab.savefig(path + 'residuals_fitfunction.png')

lab.clf()
fig = plt.figure()

c = np.array([])
for k in range(NUMCHIPS):
  c = np.append(c,np.fromfile(path + "chip_%i.csv" %(k+1), sep="\t"))
d = c.reshape((-1,14))

#Plotting used objects
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(d[:,12],d[:,13],'k.')
appearance('','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'nocolorbar', 'camgrid', ax1, 'sub')

#Plotting applied contour fit
ax2 = fig.add_subplot(2, 2, 2)
xi = lab.linspace(min(b[:,0]), max(b[:,0]))
yi = lab.linspace(min(b[:,1]), max(b[:,1]))
zi = lab.griddata(b[:,0], b[:,1], b[:,2], xi, yi)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.contourf(xi, yi, zi, 10)
axbar = ax2.contourf(xi, yi, zi, 10)
appearance('','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax2, 'sub')


ax3 = fig.add_subplot(2, 2, 3)
zi = lab.griddata(b[:,0], b[:,1], b[:,3], xi, yi)
plt.setp(ax3.get_xticklabels(), visible=True)
plt.setp(ax3.get_yticklabels(), visible=True)
ax3.contourf(xi, yi, zi, 10)
axbar = ax3.contourf(xi, yi, zi, 10)
appearance('Xpos','Ypos','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax3, 'sub')

ax4 = fig.add_subplot(2, 2, 4)
zi = lab.griddata(b[:,0], b[:,1], b[:,4], xi, yi)
plt.setp(ax4.get_xticklabels(), visible=True)
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.contourf(xi, yi, zi, 10)
axbar = ax4.contourf(xi, yi, zi, 10)
appearance('Xpos','','', [CAMXMIN, CAMXMAX], [CAMYMIN, CAMYMAX], 'nogrid', 'colorbar', 'camgrid', ax4, 'sub')
lab.savefig(path + 'camera.png')

os.popen("echo Ende plotten: `date` >> ~/dauer.txt")