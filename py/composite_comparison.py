#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"



from methods import latexify, format_axes, gauss, hist

import numpy as np
import matplotlib.pylab as pl
import seaborn as sns; sns.set_style('ticks')
#cmap = sns.cubehelix_palette(n_colors=6, start=1, rot=0.2, gamma=1.0, hue=0.8, light=0.85, dark=0.15, reverse=True, as_cmap=False)
cmap = sns.color_palette("cubehelix", 6)


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


if __name__ == '__main__':
    ratio = (1.0 + np.sqrt(5.0))/2.0
    latexify(columns=2)
    fig = pl.figure()
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    ax = fig.add_subplot(111)
    ax2 = pl.axes([0.16, 0.22, .4, .4], axisbg='white')

    #Selsing et al. 2015
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Selsing2015.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm1 = np.median(flux[norm_reg1])

    norm_reg2 = (wave > 3850) & (wave < 3875)
    norm2 = np.median(flux[norm_reg2])

    from scipy.interpolate import InterpolatedUnivariateSpline as spline
    f = spline(wave[np.where(flux != 0)], flux[np.where(flux != 0)])
    flux = f(wave)


    nbins = len(wave)
    log_binned_wl = np.array(hist(wave,[min(wave),max(wave)], int(2*nbins),'log'))
    from scipy.interpolate import InterpolatedUnivariateSpline
    sps = InterpolatedUnivariateSpline(wave, flux)
    flux = medfilt(sps(log_binned_wl) , 9)
    sps = InterpolatedUnivariateSpline(wave, err)
    err = medfilt(sps(log_binned_wl) , 9)
    wave = log_binned_wl
    wave = wave[np.where(wave <= 11400)]
    flux = flux[np.where(wave <= 11400)]
    err = err[np.where(wave <= 14000)]
    positive = (flux - err > 0 )
    ax.plot(wave[1::1], flux[1::1], label='This work', zorder=5, lw = 0.75, color = cmap[0], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[0], rasterized=True)
    ax2.plot(wave[1::1], flux[1::1], label = 'This work', zorder=5, lw = 0.75, color = cmap[0], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[0], rasterized=True)


    #Lusso et al. 2015
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Lusso2015.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)
    err *= (norm1 / norm)

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Lusso+15', zorder=4, lw = 0.5, alpha = 1.0, color = cmap[1], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[1], rasterized=True)
    ax2.plot(wave, flux, label = 'Lusso+15', zorder=4, lw = 0.5, alpha = 1.0, color = cmap[1], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[1], rasterized=True)




    #Glikman et al. 2006
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Glikman2006.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg2 = (wave > 3850) & (wave < 3875)
    norm = np.median(flux[norm_reg2])
    flux *= (norm2 / norm)
    err *= (norm2 / norm)

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Glikman+06', zorder=2, lw = 0.5, alpha = 1.0, color = cmap[2], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[2], rasterized=True)




    #Vanden Berk et al. 2001
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/VandenBerk2001.dat'
    wave, flux, err = read_text(filename=filename, err=True)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)
    err *= (norm1 / norm)


    flux = flux[(wave > 900)]
    err = err[(wave > 900)]
    wave = wave[(wave > 900)]

    positive = (flux - err > 0 )
    ax.plot(wave, flux, label = 'Vanden Berk+01', zorder=3, lw = 0.5, alpha = 1.0, color = cmap[3], linestyle='steps-mid')
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[3], rasterized=True)
    ax2.plot(wave, flux, label = 'Vanden Berk+01', zorder=3, lw = 0.5, alpha = 1.0, color = cmap[3], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[3], rasterized=True)





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
    ax.fill_between(wave, flux - err, flux +  err, alpha=0.2, label=r'1 $\sigma$ confidence interval', color = cmap[4], rasterized=True)
    ax2.plot(wave, flux, label = 'Telfer+02', zorder=1, lw = 0.5, alpha = 1.0, color = cmap[4], linestyle='steps-mid')
    ax2.fill_between(wave, flux - err, flux +  err, alpha=0.2, where=positive, label=r'1 $\sigma$ confidence interval', color = cmap[4], rasterized=True)




    #Francis et al. 1991
    filename = '/Users/jselsing/Work/Projects/QuasarComposite/py/data/templates/Francis1991orig.dat'
    wave, flux = read_text(filename=filename, err=False)
    norm_reg1 = (wave > 1425) & (wave < 1450)
    norm = np.median(flux[norm_reg1])
    flux *= (norm1 / norm)


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
    ax.set_xticks([100, 200, 500, 1000, 2000, 5000, 10000])
    ax.get_xaxis().tick_bottom()
    # ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50])
    ax.set_yticks([0.3, 1, 3, 10, 30, 100, 300])




    ax.set_xlim((300, 13000))
    ax.set_ylim((0.15, 110))

    ax2.set_xlim((700, 1400))
    ax2.set_ylim((3*3, 9*3))


    # set the linewidth of each legend object
    leg = ax.legend()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)

    format_axes(ax)
    format_axes(ax2)
    for axis in [ax2.xaxis, ax2.yaxis]:
        axis.set_tick_params(direction='in')


    pl.tight_layout()
    pl.savefig('../documents/figs/composite_comparison.pdf', rasterized=True)
    pl.show()