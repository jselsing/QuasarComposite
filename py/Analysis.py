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


#
#

def main():
    from matplotlib import rc_file
    rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
    import numpy as np
    import matplotlib.pyplot as pl



    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'

    data_file = np.genfromtxt(root_dir+'Composite.dat')


    wl = data_file[:,0]
    mean = data_file[:,1]
    err_mean = data_file[:,2]
    wmean = data_file[:,3]
    err_wmean = data_file[:,4]
    geo_mean = data_file[:,5]
    median = data_file[:,6]

    # Smoothing for presentation purposes
    from gen_methods import smooth,medfilt
    from xshoo.binning import binning1d
    # bins = 5
    # wl = binning1d(wl, bins)
    # mean, std = binning1d(mean, bins, err=errofmean)
    # mean = medfilt(mean, 5)
    # errofmean = medfilt(errofmean, 5)
    # wmean = medfilt(wmean, 5)
    # errofwmean = medfilt(errofwmean, 5)
    # geo_mean = medfilt(geo_mean, 5)
    # median = medfilt(median, 5)
    # median = smooth(np.array(median), window_len=5, window='hanning')
    # wmean = smooth(np.array(wmean), window_len=5, window='hanning')



    #Fitting power laws
    from scipy import optimize

    def power_law(x_tmp, a_tmp, k_tmp):

        tmp = a_tmp * x_tmp ** k_tmp
        return tmp

    def power_law2(x_tmp, a1_tmp, k1_tmp, a2_tmp, k2_tmp):

        tmp1 = power_law(x_tmp, a1_tmp,k1_tmp)[x_tmp<5000]

        tmp2 = power_law(x_tmp, a2_tmp,k2_tmp)[x_tmp>= 5000]

        return np.concatenate((tmp1,tmp2))

    par_guess = [1, -1.54, 1, -1.46]


    # pl.plot(wl, power_law2(wl, par_guess))
    # pl.show()
    wmean[np.where(np.isnan(wmean) == True)] = 0

    mask = (wl > 1350) & (wl < 1365) | (wl > 4200) & (wl < 4230) | (wl > 6005) & (wl < 6035) | (wl > 7160) & (wl < 7180)
    print(mask)
    popt, pcov = optimize.curve_fit(power_law2, wl[mask], wmean[mask], p0=par_guess)
    print(*popt)

    #Plotting
    pl.plot(wl, power_law2(wl, *popt))
    # pl.plot(wl, wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # pl.plot(wl, mean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter mean composite')
    # pl.plot(wl, geo_mean, color = 'red', lw = 0.5, linestyle = 'steps-mid', label='X-shooter geometric composite')
    pl.plot(wl, median, color = 'green', lw = 0.5, linestyle = 'steps-mid', label='X-shooter median composite')
    # pl.plot(wl, err_mean, color = 'black', lw = 0.5, linestyle = 'steps-mid')

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


    # ------------------------- Auxilliary products -------------------------
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


    pl.ylim((ymin,ymax))
    pl.xlim((0.95*min(wl), 1.01*max(wl)))
    pl.xlabel(r'Rest Wavelength  [\AA]')
    pl.ylabel(r'Normalised FLux [unitless]')
    # pl.semilogy()
    # pl.semilogx()
    pl.loglog()
    pl.legend()
    # pl.tight_layout()
    fig = pl.gcf()
    fig.set_size_inches(10,3)
    fig.savefig("../documents/Composite.pdf", dpi= 150)
    pl.show()



if __name__ == '__main__':
    main()