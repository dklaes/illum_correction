#!/usr/bin/env python

# Python module for photometric calibration.
# It needs the Python modules nmpfit, numpy and matplotlib 
# to be installed.

# 29.06.2013 Plot markers are adjusted according to the number of data
#            points (the old crosses were too big to see features if a
#            lot of data points were present).
# 05.06.2013 For the fit rejection trimmed means andd sigmas are estimated
#            now. We saw that outliers could dominate the statistics.
#            There is currently a hardcoded value that 10% of the objects
#            are rejected on both ends of the distribution.
# 26.11.2012 The night argument (command line) was never used. I removed it.
#
# 27.06.2012 Plotting did not work properly if the photometric solution
#            did not provide meaningful error estimates (Null value).
#            This probably originates from numerical instabilities and
#            should only occur in exceptional circumstances. I fixed it
#            with a workaround at the moment.
# 26.01.2012 The script did not work if an interactive window display
#            was not available - also in pure batch processing - fixed.
# 10.02.2010 Print results (ZPs, ext coeff. etc.) with two valid digits only
# 18.05.2009 Remove filtering for points outside plotting
#            range. matplotlib's axis() does that just fine.
# 15.05.2009 Use numpy, use matplotlib instead of ppgplot, introduced
#            batch mode
# 03.03.2005 Fixed a serious bug in the rejection loop. Instead
#            of using the remaining points we always used all points
#            and rejected points until the original fit matched the data
# 15.02.2005 Fixed the range of the y-axes in the plots to more
#            sensible values
# 14.02.2005 Fixed a bug when more paramters were fitted than
#            data points were present
#            We now rescale points to the airmass/color at which
#            they are plotted (zero)
#            Check that label is set
# 10.12.2004 Now takes a new argument "label" to be
#            used as axis label in the color plot

import copy
import getopt
import string
import sys

import numpy as np
import matplotlib
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.mathtext as mathtext
#import matplotlib.pyplot as plt
import matplotlib.artist as artist
import matplotlib.image as image
from matplotlib.widgets import RadioButtons
# note that we import matplotlib.plplot later in the
# main program. We can do this only if we know about
# batchmode processing.

import nmpfit

# Default: no acceptable solution
selVal = 0 

    
class ItemProperties:
    def __init__(self, fontsize=14, labelcolor='black', bgcolor='yellow',
                 alpha=1.0):
        self.fontsize = fontsize
        self.labelcolor = labelcolor
        self.bgcolor = bgcolor
        self.alpha = alpha

        self.labelcolor_rgb = colors.colorConverter.to_rgb(labelcolor)
        self.bgcolor_rgb = colors.colorConverter.to_rgb(bgcolor)


class MenuItem(artist.Artist):
    parser = mathtext.MathTextParser("Bitmap")
    padx = 5
    pady =5
    def __init__(self, fig, labelstr, props=None, hoverprops=None,
                 on_select=None):
        artist.Artist.__init__(self)

        self.set_figure(fig)
        self.labelstr = labelstr

        if props is None:
            props = ItemProperties()

        if hoverprops is None:
            hoverprops = ItemProperties()

        self.props = props
        self.hoverprops = hoverprops


        self.on_select = on_select

        x, self.depth = self.parser.to_mask(
            labelstr, fontsize=props.fontsize, dpi=fig.dpi)

        if props.fontsize!=hoverprops.fontsize:
            raise NotImplementedError('support for different font sizes not implemented')


        self.labelwidth = x.shape[1]
        self.labelheight = x.shape[0]

        self.labelArray = np.zeros((x.shape[0], x.shape[1], 4))
        self.labelArray[:,:,-1] = x/255.

        self.label = image.FigureImage(fig, origin='upper')
        self.label.set_array(self.labelArray)

        # we'll update these later
        self.rect = patches.Rectangle((0,0), 1,1)

        self.set_hover_props(False)

        fig.canvas.mpl_connect('button_release_event', self.check_select)

    def check_select(self, event):
        over, junk = self.rect.contains(event)
        if not over:
            return

        if self.on_select is not None:
            self.on_select(self)

    def set_extent(self, x, y, w, h):
        self.rect.set_x(x)
        self.rect.set_y(y)
        self.rect.set_width(w)
        self.rect.set_height(h)

        self.label.ox = x+self.padx
        self.label.oy = y-self.depth+self.pady/2.

        self.rect._update_patch_transform()
        self.hover = False

    def draw(self, renderer):
        self.rect.draw(renderer)
        self.label.draw(renderer)

    def set_hover_props(self, b):
        if b:
            props = self.hoverprops
        else:
            props = self.props

        r, g, b = props.labelcolor_rgb
        self.labelArray[:,:,0] = r
        self.labelArray[:,:,1] = g
        self.labelArray[:,:,2] = b
        self.label.set_array(self.labelArray)
        self.rect.set(facecolor=props.bgcolor, alpha=props.alpha)

    def set_hover(self, event):
        'check the hover status of event and return true if status is changed'
        b,junk = self.rect.contains(event)

        changed = (b != self.hover)

        if changed:
            self.set_hover_props(b)

        self.hover = b
        return changed


class Menu:
    def __init__(self, fig, menuitems):
        self.figure = fig
        fig.suppressComposite = True

        self.menuitems = menuitems
        self.numitems = len(menuitems)

        maxw = max([item.labelwidth for item in menuitems])
        maxh = max([item.labelheight for item in menuitems])

        totalh = self.numitems*maxh + (self.numitems+1)*2*MenuItem.pady
        
        x0 = 350
        y0 = 75

        width = maxw + 2*MenuItem.padx
        height = maxh+MenuItem.pady

        for item in menuitems:
            left = x0
            bottom = y0-maxh-MenuItem.pady
            item.set_extent(left, bottom, width, height)
            fig.artists.append(item)
            y0 -= maxh + MenuItem.pady
        fig.canvas.mpl_connect('motion_notify_event', self.on_move)

    def on_move(self, event):
        draw = False
        for item in self.menuitems:
            draw = item.set_hover(event)
            if draw:
                self.figure.canvas.draw()
                break
            

# estimate trimmed mean and standard deviations from
# an array:
def trimmed_mean_std(data, lowerpercentile, upperpercentile):
    data = np.array(data)
    data.sort()
    low = int(lowerpercentile * len(data))
    high = int((1. - upperpercentile) * len(data))
    return (data[low:high].mean(), data[low:high].std(ddof=0))
        
def phot_funct(p, fjac=None, airmass=None, color=None, y=None, err=None):
    model = p[0] + p[1]*airmass + p[2]*color

    status = 0
    return([status, model-y])


def readInput(file):
    f = open(file, "r")

    instMagList = []
    stdMagList = []
    colList = []
    airmassList = []

    for line in f.readlines():
        instMag, stdMag, col, airmass = string.split(line)
        instMagList.append(float(instMag))
        stdMagList.append(float(stdMag))
        colList.append(float(col))
        airmassList.append(float(airmass))
    f.close()

    instMag = np.array(instMagList)
    stdMag = np.array(stdMagList)
    data = stdMag - instMag
    airmass = np.array(airmassList)
    color = np.array(colList)

    return data, airmass, color


def photCalib(data_save, airmass_save, color_save, p, maxSigIter=50):

    save_len = len(data_save)
    
    parinfo = [{"value": p[0], "fixed": 0},
               # The extinction coefficient has to be negative
               {"value": p[1], "fixed": 0, "limited": [0,1],
                "limits": [-99, 0]},
               {"value": p[2], "fixed": 0}]
    

    solutions = []

    for fixedList in [[], [1], [1,2]]:
        airmass = copy.copy(airmass_save)
        color = copy.copy(color_save)
        data = copy.copy(data_save)

        fa = {"airmass": airmass,
              "color": color,
              "y": data}

        for j in range(len(parinfo)):
            if j in fixedList:
                print "Element", j, "is fixed at", p[j]
                parinfo[j]["fixed"] = 1
            else:
                parinfo[j]["fixed"] = 0

        for i in range(maxSigIter):
            old_len = len(data)
            fa = {"airmass": airmass,
                  "color": color,
                  "y": data}
            m =  nmpfit.mpfit(phot_funct, functkw=fa,
                              parinfo=parinfo,
                              maxiter=1000, quiet=1)
            
            if (m.status <= 0):
                print 'error message = ', m.errmsg
                condition = np.zeros(len(data))
                break

            airmass = copy.copy(airmass_save)
            color = copy.copy(color_save)
            data = copy.copy(data_save)

            # Compute a 3 sigma rejection criterion
            condition = threeSigmaCond(m.params, data_save, data,
                                       airmass_save, airmass,
                                       color_save, color)
            
            # Keep everything (from the full data set!) that is within
            # the 3 sigma criterion
            data = np.compress(condition, data_save)
            airmass = np.compress(condition, airmass_save)
            color = np.compress(condition, color_save)
            new_len = len(data)
            
            if float(new_len)/float(save_len) < 0.5:
                print "Rejected more than 50% of all measurements."
                print "Aborting this fit."
                break
            
            # No change
            if new_len == old_len:
                print "Converged! (%d iterations)" % (i+1, )
                print "Kept %d/%d stars." % (new_len, save_len)
                break

        solutions.append([m.params, m.perror, condition])
        print 
    return solutions


def threeSigmaCond(p, data_save, data, airmass_save, airmass, color_save,
                   color):
    if len(data_save) > 1:
        mo = p[0] + p[1]*airmass + p[2]*color
        mo_save = p[0] + p[1]*airmass_save + p[2]*color_save
        dm = data-mo
        dm_save = data_save - mo_save
        #mean =  np.sum(dm)/len(dm)
        #mmsigma = np.sqrt(np.sum(np.power(mean-dm, 2)) / (len(dm) - 1))

        # estimate trimmed mean and sigma values. We trim 10% of the
        # data in the upper and lower parts
        mean, sigma = trimmed_mean_std(dm, 0.1, 0.1)

        condition = np.less(np.fabs(dm_save), 3. * sigma)
    else:
        condition = np.zeros(len(data_save))
    
    return condition


def makePlots(data, airmass, color, outfile, solutions, label):
    global selVal
    file = outfile+".png"

    fig = plt.figure(num=None, figsize=(15, 13))
    plt.subplots_adjust(top=0.97, right=0.99, bottom=0.15, wspace=0.15)

    for i in range(3):
        result = solutions[i]

        # Airmass plot
        plt.subplot(3, 2, 2 * i + 1)
        airMin = 1
        airMax = np.sort(airmass)[-1] * 1.1
        dataAirMax = result[0][0] + result[0][1] + 1
        dataAirMin = result[0][0] + result[0][1] - 1
        dataColMax = result[0][0] + 1
        dataColMin = result[0][0] - 1
        colMinVal = np.sort(color)[0]
        if colMinVal < 0:
            colMin = colMinVal * 1.1
        else:
            colMin = colMinVal * 0.95
        colMax = np.sort(color)[-1] * 1.1
        
        if result[0] != None and result[1] != None and \
                result[0].any() and result[1].any():
            eqStr = "%d parameter fit: Mag$-$Mag(Inst) = $%.2f \pm %.2f + \
                    (%.2f \pm %.2f)$ airmass$ + "\
                    "(%.2f \pm %.2f)$ color" % \
                    (3-i, result[0][0], result[1][0], \
                     result[0][1], result[1][1], result[0][2], result[1][2])
        else:
            eqStr = "%d parameter fit not possible" % (3-i, )
        

        condition = result[2]
        goodAirmass = np.compress(condition, airmass)
        goodData = np.compress(condition, data)
        goodColor = np.compress(condition, color)
        badAirmass = np.compress(np.logical_not(condition), airmass)
        badData = np.compress(np.logical_not(condition), data)
        badColor = np.compress(np.logical_not(condition), color)
        
        # adjust plot markers and marker sizes according to the
        # number of data points:
        if (len(goodData) + len(badData)) < 100:
            markersymbol="x"
            markersize=3
        else:
            markersymbol="."
            markersize=1

        if len(goodData):
            # Rescale to zero color
            plotData = goodData - result[0][2] * goodColor
            plt.scatter(goodAirmass, plotData, color="g",
                        marker=markersymbol, s=markersize)
        if len(badData):
            plotData = badData - result[0][2] * badColor
            plt.scatter(badAirmass, plotData, color="r",
                        marker=markersymbol, s=markersize)


        a = np.arange(1, airMax, 0.01)
        m = result[0][0] + result[0][1] * a
        plt.plot(a, m, linestyle='solid')
        fixenv([1, airMax], [dataAirMin, dataAirMax],
               eqStr, label=["Airmass", "Mag - Mag(Inst)"])
        
        # Color Plot
        plt.subplot(3, 2, 2 * i + 2)

        if len(goodData):
            # Rescale to zero airmass 
            plotData = goodData-result[0][1]*goodAirmass
            plt.scatter(goodColor, plotData, color="g",
                        marker=markersymbol, s=markersize)
        if len(badData):
            plotData = badData-result[0][1]*badAirmass
            plt.scatter(badColor, plotData, color="r",
                        marker=markersymbol, s=markersize)

        a = np.arange(colMin, colMax, 0.01)
        m = result[0][0] + result[0][2] * a
        plt.plot(a, m, linestyle='solid')
        fixenv([colMin, colMax] ,
               [dataColMin, dataColMax],
               eqStr, label=[label, "Mag - Mag(Inst)"])

    plt.savefig(file)

    axcolor = 'lightgoldenrodyellow'
    rax = plt.axes([0.1, 0, 0.15, 0.12], axisbg=axcolor)
    radio = RadioButtons(rax, ("No acceptable solution",
                               "3 parameter fit",
                               "2 parameter fit",
                               "1 parameter fit"))
    radio.on_clicked(updateSel)
    props = ItemProperties(labelcolor='black', bgcolor='lightgoldenrodyellow',
                           fontsize=15, alpha=1)
    hoverprops = ItemProperties(labelcolor='white', bgcolor='red',
                                fontsize=15, alpha=1)

    menuitems = []
    for label in ('Done',):
        item = MenuItem(fig, label, props=props, hoverprops=hoverprops,
                        on_select=on_done)
        menuitems.append(item)
    menu = Menu(fig, menuitems)
    return

def on_done(item):
    sys.exit(selVal)


def updateSel(label):
    global selVal
    selDict = {"3 parameter fit": 1, "2 parameter fit": 2,
               "1 parameter fit": 3, "No acceptable solution": 0}
    selVal = selDict[label]
    plt.draw()
    return
    

def fixenv (xrange=[0,1], yrange=[0,1], fname="none", ci = 1,
            label=["x", "y"]):                              
    plt.xlabel(label[0], size="small")
    plt.ylabel(label[1], size="small")
    plt.title(fname, size="small")
    plt.axis([xrange[0], xrange[1], yrange[0], yrange[1]]) # set axis ranges.
    return


def saveResults(file, solutions):
    f = open(file+".asc", "w")
    for result in solutions:
        if np.sometrue(result[2]):
            f.write("%.2f %.2f %.2f %.2f %.2f %.2f\n" % \
                        (result[0][0], result[0][1], result[0][2], result[1][0], result[1][1], result[1][2]))
        else:
            f.write("-1 -1 -1 -1 -1 -1\n")
    f.close
    return
                


def usage():
    print "Usage:"
    print "photo_abs.py [-b] -i input -f filter -n GABODSID - e ext. coeff. -c color coeff. -o output -l label"
    print
    print "    -b, --batch               Run in batch mode, no interactive selection"
    print "    -i, --input=STRING        Input file, must have 4 columns: Instrumental Mag, Standard Mag, Color, Airmass"
    print "    -o, --output=STRING       Output file basename"
    print "    -e, --extinction=FLOAT    Default value of extinction coefficient for one/two parameter fit"
    print "    -c, --color=FLOAT         Default value of color term for one parameter fit"
    print "    -l, --label=STRING        Label for color axis (e.g. B-V)"
    print
    print "Author:"
    print "    Joerg Dietrich <jdietric@eso.org>"
    print
    print "Modified by Dominik Klaes <dklaes@astro.uni-bonn.de>: smaller dots in plots and errors from the fit are also"
    print "written into the result file."
    print
    return
    


if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "bi:n:o:e:c:l:",
                                   ["batch", "input=", "night=", "extinction=",
                                    "color=", "output=", "label="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    infile = night = extcoeff = colcoeff = outfile = label = batchmode = None
    for o, a in opts:
        if o in ("-b", "--batch"):
            batchmode = True
        elif o in ("-i", "--input"):
            infile = a
        elif o in ("-o", "--output"):
            outfile = a
        elif o in ("-e", "--extinction"):
            extcoeff  = float(a)
        elif o in ("-c", "--color"):
            colcoeff = float(a)
        elif o in ("-l", "--label"):
            label = a
        else:
            print "option:", o
            usage()
            sys.exit(2)
            
    if not infile or not outfile or \
           extcoeff==None or colcoeff==None or label==None:
        print infile, night, outfile, coeff, color
        usage()
        sys.exit(2)

    # we import pyplot.pyplot here after we know whether we
    # need batch processing. If we run in batch mode we disable
    # use of the interactive window display.
    if batchmode:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    data, airmass, color = readInput(infile)
    solutions = photCalib(data, airmass, color, [24, extcoeff, colcoeff])
    makePlots(data, airmass, color, outfile, solutions, label)
    saveResults(outfile, solutions)

    if not batchmode:
        plt.show()
    sys.exit(selVal)
