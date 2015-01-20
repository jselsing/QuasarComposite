#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SYNOPSIS

    TODO helloworld [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    TODO: Name <name@example.org>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

VERSION

    $
"""

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')


import numpy as np
import glob
from scipy import interpolate
import matplotlib.pylab as pl
from unred import ccm_unred



if __name__ == '__main__':
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/Telluric_corrected_science.dat')
    redshifts = [1.1250, 1.9798, 1.5826, 1.8458, 1.5123, 2.0998, 1.3092]
    ebv = [0.0253, 0.0201, 0.0279, 0.0268, 0.0330, 0.0274, 0.0378]
    sdss_data_files = np.array([np.genfromtxt(i) for i in sdssobjects])

    wl = np.array([sdss_data_files[i][:,0] / (1 + redshifts[i]) for i in range(len(sdssobjects))])

    flux = np.array([sdss_data_files[i][:,1] for i in range(len(sdssobjects))])
    fluxerr = np.array([sdss_data_files[i][:,2] for i in range(len(sdssobjects))])


    #de-reddening
    # wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    # flux = np.array([ccm_unred(wl_obs[0],sdss_data_files[i][:,1][0], ebv[i]) for i in range(len(sdssobjects))])
    # fluxerr = np.array([ccm_unred(wl_obs[0],sdss_data_files[i][:,2][0], ebv[i]) for i in range(len(sdssobjects))])


    # Interpolate to a common wavelength:
    short = []
    tall = []
    for i in wl:
        short.append(min(i))
        tall.append(max(i))
    short = min(short)
    tall = max(tall)

    step = 0.2 #CDELT
    #wl_new = (np.arange(((tall-short)/step))*step+short)
    wl_new = np.arange(short, tall, step)
    print(short, tall, len(wl_new))
    flux_new = np.zeros((len(redshifts),len(wl_new)))
    fluxerr_new = np.zeros((len(redshifts),len(wl_new)))
    for n in range(np.shape(wl)[0]):
        f = interpolate.interp1d(wl[n],flux[n],kind='linear',bounds_error = False, fill_value=0.)
        g = interpolate.interp1d(wl[n],fluxerr[n],kind='linear',bounds_error = False, fill_value=1.)
        mask = (wl_new > 7000) & (wl_new < 7500)
        norm = np.median(f(wl_new)[mask])
        flux_new[n] = f(wl_new)/norm
        fluxerr_new[n] = g(wl_new)/norm
    #     pl.plot(wl_new, flux_new[n], lw=0.2)
    # pl.show()
    wl = wl_new
    flux = flux_new
    fluxerr = fluxerr_new



    #------------------------- Combination -------------------------
    # Weighted average:
    weight = 1./(np.array(fluxerr)**2)
    # mean = np.average(flux, axis = 0, weights = weight)
    mean = np.median(flux, axis = 0)
    import scipy
    import scikits.bootstrap as bootstrap
    CI_low = np.zeros(np.shape(mean))
    CI_high = np.zeros(np.shape(mean))
    for i,j in enumerate(flux.transpose()[:-1]):
        CI_low[i], CI_high[i] = bootstrap.ci(data=j, statfunction=np.average, n_samples=100, method='bca', epsilon=0.001)
        if (i % 1000==0):
            print(i)
    boot_mean = (CI_high + CI_low)/2
    # from gen_methods import mytotal
    # mean, errofmean, lo, hi = mytotal(flux, axis = 2, type='median')
    from xshoo.binning import binning1d
    wl_ci = binning1d(wl, 5)
    CI_low = binning1d(CI_low, 5)
    CI_high = binning1d(CI_high, 5)

    pl.plot(wl, mean)
    pl.plot(wl, boot_mean, color='black', lw=0.3)
    pl.fill_between(wl_ci, CI_low,CI_high, color= 'red', alpha = 0.2)
    pl.show()


    errofmean = np.sqrt(1./np.sum(np.array(fluxerr)**-2.,axis=0))


    fluxerr_new = []
    for j , (k, l) in enumerate(zip(mean[:-1],errofmean[:-1])):
        # if k > 3 * mean[j-1] and k > 0:
        #    fluxerr_new.append(abs(100*errofmean[j+1]))
        if k < 0.75 * mean[j-1] and k > 0:
            fluxerr_new.append(abs(100*errofmean[j+1]))
        else:
            fluxerr_new.append(abs(errofmean[j+1]))
    fluxerr_new.append(0)
    from gen_methods import smooth
    errofmean = smooth(np.array(fluxerr_new), window_len=2, window='hanning')





    # from xshoo.binning import binning1d
    # bins = 5
    # mean, std = binning1d(mean, bins, err=errofmean)
    # wl = binning1d(wl, bins)
    # mean = smooth(np.array(mean), window_len=5, window='hanning')

    pl.plot(wl, mean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter composite')
    #pl.plot(wl_new, errofmean, color = 'black', lw = 0.5, linestyle = 'steps-mid')

    #Overplot lines
    fit_line_positions = np.genfromtxt('linelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)

    for p in range(len(fit_line_positions)):
        xcoord = linelist[p]
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(mean[mask])
        pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
        pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.1 ),fontsize='x-small')

    #------------------------- Auxilliary products -------------------------
    # Johans Composite:
    vandenberk = np.loadtxt('Glikman.data')
    mask = (wl > 3600) & (wl < 3700)
    maskVB = (vandenberk[:,0] > 3600) & (vandenberk[:,0] < 3700)
    norm = np.median(vandenberk[:,1][maskVB]) / np.median(mean[mask])
    # Vanden Berk:
    vandenberk2 = np.loadtxt('composite.txt')
    mask = (wl > 3600) & (wl < 3700)
    maskVB = (vandenberk2[:,0] > 3600) & (vandenberk2[:,0] < 3700)
    norm2 = np.median(vandenberk2[:,1][maskVB]) / np.median(mean[mask])
    pl.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Glikman Composite')
    pl.plot(vandenberk2[:,0],vandenberk2[:,1]/norm2,drawstyle='steps-mid',label='Vanden Berk Composite')

    pl.xlabel(r'Rest Wavelength  [\AA]')
    pl.ylabel(r'Normalised FLux [unitless]')
    pl.semilogy()
    pl.legend()
    pl.show()