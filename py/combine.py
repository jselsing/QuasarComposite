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
from unred import ccm_unred,cardelli_reddening


def main():
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/Telluric_corrected_science.dat')
    redshifts = [1.1250, 1.9798, 1.5826, 1.8458, 1.5123, 2.0998, 1.3092]
    ebv = [0.0253, 0.0201, 0.0279, 0.0268, 0.0330, 0.0274, 0.0378]
    sdss_data_files = np.array([np.genfromtxt(i) for i in sdssobjects])

    wl = np.array([sdss_data_files[i][:,0] / (1 + redshifts[i]) for i in range(len(sdssobjects))])
    wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    flux = np.array([sdss_data_files[i][:,1] for i in range(len(sdssobjects))])
    fluxerr = np.array([sdss_data_files[i][:,2] for i in range(len(sdssobjects))])
    bp_map = np.array([sdss_data_files[i][:,3] for i in range(len(sdssobjects))])

    # Interpolate to a common wavelength:
    short = []
    tall = []
    for i in wl:
        short.append(min(i))
        tall.append(max(i))
    short = min(short)
    tall = max(tall)

    step = 0.2 #CDELT
    wl_new = np.arange(short, tall, step)
    print(short, tall, len(wl_new))
    flux_new = np.zeros((len(redshifts),len(wl_new)))
    fluxerr_new = np.zeros((len(redshifts),len(wl_new)))
    bp_map_new = np.zeros((len(redshifts),len(wl_new)))
    for n in range(np.shape(wl)[0]):
        #de-reddening
        flux[n] = cardelli_reddening(wl_obs[n], flux[n], ebv[n])
        fluxerr[n] = cardelli_reddening(wl_obs[n], fluxerr[n], ebv[n])
        #Interpolate
        f = interpolate.interp1d(wl[n],flux[n],kind='linear',bounds_error = False, fill_value=0.)
        g = interpolate.interp1d(wl[n],fluxerr[n],kind='linear',bounds_error = False, fill_value=1.)
        h = interpolate.interp1d(wl[n],bp_map[n],kind='linear',bounds_error = False, fill_value=1.)
        mask = (wl_new > 7000) & (wl_new < 7500)
        norm = np.median(f(wl_new)[mask])
        flux_new[n] = f(wl_new)/norm
        fluxerr_new[n] = g(wl_new)/norm
        bp_map_new[n] = h(wl_new)
    #     pl.plot(wl_new, flux_new[n], lw=0.2)
    # pl.show()
    wl = wl_new
    flux = flux_new
    fluxerr = fluxerr_new
    bp_map = bp_map_new


    #------------------------- Combination -------------------------
    # Weighted average:
    wmean = np.zeros(np.shape(wl))
    mean = np.zeros(np.shape(wl))
    geo_mean = np.zeros(np.shape(wl))
    median = np.zeros(np.shape(wl))
    errofmean = np.zeros(np.shape(wl))
    import scipy
    import scikits.bootstrap as bootstrap
    CI_low = np.zeros(np.shape(mean))
    CI_high = np.zeros(np.shape(mean))
    for i, k in enumerate(flux.transpose()):
        mask = np.where(bp_map.transpose()[i] == 0)
        if len(k[mask]) != 0:
            #Weighted mean
            weight = 1./(np.array(fluxerr.transpose()[i][mask])**2)
            wmean[i] = np.average(k[mask], axis = 0, weights = weight)
            errofmean[i] = np.sqrt(1./np.sum(np.array(fluxerr.transpose()[i][mask])**-2.,axis=0))


            #Mean
            mean[i] = np.mean(k[mask])
            # errofmean[i] = np.std(fluxerr.transpose()[i][mask])


            #Geometric mean
            from scipy.stats.mstats import gmean
            geo_mean[i] = gmean(k[mask])

            #Median
            median[i] = np.median(k[mask])
            # Bootstrapping to get confidence intervals
            # if len(k[mask]) > 1:
            #     CI_low[i], CI_high[i] = bootstrap.ci(data=k[mask], statfunction=np.median, n_samples=100, method='bca')
            # else:
            #     CI_low[i], CI_high[i] = np.median(k[mask]), np.median(k[mask])
            # if (i % 1000==0):
            #     print(i)
        else:
            mean[i] = 0
            geo_mean[i] = 0
            median[i] = 0
            CI_low[i], CI_high[i] = 0, 0

    pl.plot(wl, mean/errofmean, lw=0.5, color = 'black')
    pl.xlabel(r'Rest wavelength [\AA]')
    pl.ylabel(r'S/N')
    pl.savefig("../documents/signal_to_noise.pdf", dpi= 150)
    pl.clf()

    #Calculation of unmasked weighted mean
    # weight = 1./(np.array(fluxerr)**2)
    # mean = np.average(flux, axis = 0, weights = weight)
    # errofmean = np.sqrt(1./np.sum(np.array(fluxerr)**-2.,axis=0))
    # mean = np.median(flux, axis = 0)



    # Bootstrapping to get confidence intervals
    # import scipy
    # import scikits.bootstrap as bootstrap
    # CI_low = np.zeros(np.shape(mean))
    # CI_high = np.zeros(np.shape(mean))
    # for i,j in enumerate(flux.transpose()[:-1]):
    #     CI_low[i], CI_high[i] = bootstrap.ci(data=j, statfunction=scipy.average, n_samples=100, method='abc', epsilon=weight.transpose()[i])
    #     if (i % 1000==0):
    #         print(i)
    # boot_mean = (CI_high + CI_low)/2
    # from xshoo.binning import binning1d
    # wl_ci = binning1d(wl, 5)
    # CI_low = binning1d(CI_low, 5)
    # CI_high = binning1d(CI_high, 5)
    # pl.plot(wl, median, color='black', lw=0.3)
    # pl.plot(wl, boot_mean, color='red', lw=0.6)
    # pl.fill_between(wl_ci, CI_low , CI_high, color= 'red', alpha = 0.2)
    # pl.show()
    # pl.plot(wl, mean, color='black', lw = 0.5)
    # pl.show()



    #Calculating the number of spectra that goes into the composite and appending to n_spec.
    n_spec = np.zeros(np.shape(mean))
    spec = []
    for i, n in enumerate(flux.transpose()):
        mask = np.where(bp_map.transpose()[i] == 0)
        n_spec[i] = len((n[mask])[np.where((n[mask]) != 0)])
        if len(n[np.where(n != 0)]) == 7:
            spec.append(n / np.median(n))

    pl.plot(wl, n_spec, lw = 0.5, color = 'black')
    pl.xlabel(r'Rest wavelength [\AA]')
    pl.ylabel(r'Number of spectra')
    pl.savefig("../documents/number_spec.pdf", dpi= 150)
    pl.clf()

    # Checking for normality
    # from matplotlib import pyplot as plt
    # import matplotlib.mlab as mlab
    # import scipy.stats as stats
    # p_val = []
    # for i, k in enumerate(spec[10000:10110]):
    #     # pl.plot(n, '.' , hold=True)
    #     k = k / np.median(k)
    #     k = np.hstack((k,np.mean(k)))
    #     #print(k)
    #     p_val.append((stats.normaltest(k)[1]))
    #     n, bins, patches = plt.hist(k, 10, hold=True)
    #     mu = np.mean(k)
    #     sigma = np.std(k)
    #     plt.plot(bins, mlab.normpdf(bins, mu, sigma), hold=True)
    # print(np.mean(p_val))
    # pl.xlabel(r'Normalised value [input/median(input)]')
    # pl.ylabel(r'Arbitrary scale')
    # pl.savefig("../documents/normality.pdf", dpi= 150)



    #Smoothing for presentation purposes
    from gen_methods import smooth,medfilt
    # from xshoo.binning import binning1d
    # bins = 5
    # wl = binning1d(wl, bins)
    # mean, std = binning1d(mean, bins, err=errofmean)
    mean = medfilt(mean, 5)
    wmean = medfilt(wmean, 5)
    geo_mean = medfilt(geo_mean, 5)
    median = medfilt(median, 5)
    #mean = smooth(np.array(mean), window_len=5, window='hanning')



    #Plotting
    pl.plot(wl, wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # pl.plot(wl, mean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter mean composite')
    # pl.plot(wl, geo_mean, color = 'red', lw = 0.5, linestyle = 'steps-mid', label='X-shooter geometric composite')
    # pl.plot(wl, median, color = 'green', lw = 0.5, linestyle = 'steps-mid', label='X-shooter median composite')
    #pl.plot(wl_new, errofmean, color = 'black', lw = 0.5, linestyle = 'steps-mid')

    #Overplot lines
    fit_line_positions = np.genfromtxt('plotlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)

    ymin = 0.3
    ymax = 155

    for p in range(len(fit_line_positions)):
        xcoord = linelist[p]
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(mean[mask])
        pl.axvline(x=xcoord, ymin=(np.log(y_val) + 1.30) / (( np.log(ymax) - np.log(ymin))), color='red',linestyle='dashed', lw=0.75)
        pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 2.0 ),fontsize='x-small', rotation = 90)


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
    # pl.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Glikman Composite')
    # pl.plot(vandenberk2[:,0],vandenberk2[:,1]/norm2,drawstyle='steps-mid',label='Vanden Berk Composite')


    pl.ylim((ymin,ymax))
    pl.xlim((0.95*short, 1.01*tall))
    pl.xlabel(r'Rest Wavelength  [\AA]')
    pl.ylabel(r'Normalised FLux [unitless]')
    pl.semilogy()
    pl.legend()
    pl.tight_layout()
    fig = pl.gcf()
    fig.set_size_inches(10,3)
    fig.savefig("../documents/Composite.pdf", dpi= 150)
    pl.show()





if __name__ == '__main__':
    main()