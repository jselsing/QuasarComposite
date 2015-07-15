#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


import numpy as np
import matplotlib.pylab as pl
import seaborn as sns; sns.set_style('ticks')
#cmap = sns.cubehelix_palette(n_colors=6, start=1, rot=0.2, gamma=1.0, hue=0.8, light=0.85, dark=0.15, reverse=True, as_cmap=False)
cmap = sns.color_palette("cubehelix", 6)

import matplotlib

from math import sqrt
SPINE_COLOR = 'gray'

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'Qt4Agg',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'font.size': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              #'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)


def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(True)

    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax



def read_text(filename = 'test.dat', wave_ext = 0, flux_ext = 1, err_ext = 2, err=False):
    """Small helper to easily load texts into arrays
        Relies on numpy's genfromtxt
    """
    #Load txt file
    data = np.genfromtxt(filename)
    #Slice arrays
    wave = data[:,wave_ext]
    flux = data[:,flux_ext]
    #If error is set then return also error
    if err:
        error = data[:,err_ext]
        return wave, flux, error
    return wave, flux


from gen_methods import medfilt, smooth
def hist(rawData,xRange,nBins=10,mode='lin'):

    """histogram using linear binning of supplied data

    Input:
    rawData	-- list containing data to be binned
    xRange  -- lower(incl)/upper(excl) boundary for numerical values
    nBin    -- desired number of bins (default =10)
    mode	-- binning type (possible choices: lin, log)

    Returns: (nothing)
    """
    from math   import sqrt,floor,log,exp
    h = [0]*nBins
    xMin=float(xRange[0])
    xMax=float(xRange[1])

    if mode == 'lin':
        dx = (xMax-xMin)/nBins
        def binId(val):   return int(floor((val-xMin)/dx))
        def bdry(bin):	  return xMin+bin*dx, xMin+(bin+1)*dx
        def GErr(q,n,dx): return sqrt(q*(1-q)/(N-1))/dx

    elif mode == 'log':
        dx = log(xMax/xMin)/nBins
        def binId(val):   return int(floor(log(val/xMin)/dx))
        def bdry(bin):	  return xMin*exp(bin*dx), xMin*exp((bin+1)*dx)
        def GErr(q,n,dx): return "##"

    for value in rawData:
        if 0<=binId(value)<nBins:
          h[binId(value)] += 1

    N=sum(h)
    binned = []
    for bin in range(nBins):
        hRel   = float(h[bin])/N
        low,up = bdry(bin)
        binned.append(low)
        width  = up-low
        # print(low, up, hRel/width, GErr(hRel,N,width))
    return binned











if __name__ == '__main__':
    ratio = (1.0 + np.sqrt(5.0))/2.0
    latexify()
    fig = pl.figure(figsize=(5*ratio, 5))
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    ax = fig.add_subplot(111)
    ax2 = pl.axes([0.13, 0.2, .4, .4], axisbg='white')

    #Selsing et al. 2015
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Selsing2015.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm1 = np.median(flux[norm_reg1])

    norm_reg2 = (wave > 3850) & (wave < 3875)
    norm2 = np.median(flux[norm_reg2])


    # msk = (wave <  13000)
    # wave[msk] = smooth(wave[msk], window='hanninh', window_len=10)
    # flux[msk] = smooth(flux[msk], window='hanning', window_len=20)
    # err[msk] = err[msk]



    nbins = len(wave)
    log_binned_wl = np.array(hist(wave,[min(wave),max(wave)], int(2*nbins),'log'))
    from scipy.interpolate import InterpolatedUnivariateSpline
    sps = InterpolatedUnivariateSpline(wave, flux)
    flux = medfilt(sps(log_binned_wl) , 9)
    sps = InterpolatedUnivariateSpline(wave, err)
    err = medfilt(sps(log_binned_wl) , 9)
    wave = log_binned_wl


    positive = (flux - err > 0 )
    ax.plot(wave, flux, label='This work', zorder=5, lw = 0.75, color = cmap[0], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[0])
    ax2.plot(wave, flux, label = 'This work', zorder=5, lw = 0.75, color = cmap[0], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[0])


    #Lusso et al. 2015
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Lusso2015.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)
    err *= (norm1 / norm)

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Lusso+15', zorder=4, lw = 0.5, alpha = 1.0, color = cmap[1], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[1])
    ax2.plot(wave, flux, label = 'Lusso+15', zorder=4, lw = 0.5, alpha = 1.0, color = cmap[1], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[1])




    #Glikman et al. 2006
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Glikman2006.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg2 = (wave > 3850) & (wave < 3875)
    norm = np.median(flux[norm_reg2])
    flux *= (norm2 / norm)
    err *= (norm2 / norm)

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Glikman+06', zorder=2, lw = 0.5, alpha = 1.0, color = cmap[2], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[2])




    #Vanden Berk et al. 2001
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/VandenBerk2001.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)
    err *= (norm1 / norm)

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Vanden Berk+01', zorder=3, lw = 0.5, alpha = 1.0, color = cmap[3], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[3])
    ax2.plot(wave, flux, label = 'Vanden Berk+01', zorder=3, lw = 0.5, alpha = 1.0, color = cmap[3], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[3])





    #Telfer et al. 2002
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Telfer2002.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)
    err *= (norm1 / norm)
    mask = (wave < 2000)
    wave = wave[mask]
    flux = flux[mask]
    err = err[mask]

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Telfer+02', zorder=1, lw = 0.5, alpha = 1.0, color = cmap[4], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[4])
    ax2.plot(wave, flux, label = 'Telfer+02', zorder=1, lw = 0.5, alpha = 1.0, color = cmap[4], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[4])




    #Francis et al. 1991
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Francis1991orig.dat'
    wave, flux = read_text(filename=filename, err=False)
    #flux = flux/wave
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)

    #mask = (wave < 2000)
    #wave = wave[mask]
    #flux = flux[mask]

    ax.plot(wave, flux, label = 'Francis+91', zorder=4, lw = 1.0, alpha = 1.0, color = cmap[5], linestyle='steps-mid')
    ax2.plot(wave, flux, label = 'Francis+91', zorder=4, lw = 1.0, alpha = 1.0, color = cmap[5], linestyle='steps-mid')


    x = np.arange(0 , 15000, 0.1)
    y = x ** -(1.70)
    norm_reg1 = (x > 1425) & (x < 1450)
    norm = np.median(y[norm_reg1])
    y *= (norm1 / norm)
    ax.plot( x , y)







    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel(r'Rescaled flux density F$_\lambda$')



    ax.loglog()
    # ax.semilogy()

    # Formatting axes
    import matplotlib as mpl
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xticks([1000, 2000, 3000, 5000, 10000])
    ax.get_xaxis().tick_bottom()
    # ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50])




    ax.set_xlim((300, 13000))
    ax.set_ylim((0.1, 30))

    ax2.set_xlim((700, 1400))
    ax2.set_ylim((3, 9))
    #ax2.semilogy()


    #ax2.semilogx(data[3:8,1],data[3:8,2])
    #pl.setp(ax2, xticks=[], yticks=[])




    # set the linewidth of each legend object
    leg = ax.legend()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)

    format_axes(ax)
    format_axes(ax2)
    pl.tight_layout()
    pl.savefig('../documents/figs/composite_comparison.pdf')
    pl.show()