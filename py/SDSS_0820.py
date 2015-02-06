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
import matplotlib.pyplot as pl
from numpy.polynomial import chebyshev



if __name__ == '__main__':

    #Get spectrum
    object = []
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS0820+1306/'
    obs = np.genfromtxt(glob.glob(root_dir+'Telluric_corrected_science.dat')[0])
    wl = obs[:,0]
    flux = obs[:,1]
    fluxerr = obs[:,2]

    fluxerr_new = []
    for j, (k, l) in enumerate(zip(flux,fluxerr)):
        if k > 2 * flux[j-2] and k > 0:
            fluxerr_new.append(l + 2e-16)
        elif k < 1/2 * flux[j-2] and k > 0:
            fluxerr_new.append(l + 2e-16)
        else:
            fluxerr_new.append(l)
    from gen_methods import smooth
    fluxerr = smooth(np.array(fluxerr_new), window_len=5, window='hanning')

    sdss_wl = (obs[:,3])[np.where(obs[:,3] != 0)]
    sdss_flux = (obs[:,4])[np.where(obs[:,3] != 0)]
    redshifts = 1.1257


    # Load linelist
    fit_line_positions = np.genfromtxt('fitlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)
    linelist = linelist / (1.0 + 2.735182E-4 + 131.4182 / linelist**2 + 2.76249E8 / linelist**4)



    #Cut out fitting region
    mask = np.logical_and(wl > 10570, (wl < 10700))
    wl_fit = wl[mask]
    flux_fit = flux[mask]
    fluxerr_fit = fluxerr[mask]

    #Fit continuum and subtract
    from methods import continuum_fit
    cont, chebfit = continuum_fit(wl_fit, flux_fit, fluxerr_fit, edge_mask_len=2)
    chebfitval = chebyshev.chebval(wl, chebfit)

    #Define models to use
    from methods import voigt,gauss
    def model1(t,  amp2, sig22g, sig22l, z):
            tmp = voigt(t, abs(amp2), (1+z)*linelist[2], sig22g, sig22l)
            return tmp

    def model2(t, amp, sig2, z):
            tmp = gauss(t, abs(amp), (1+z)*linelist[2], sig2)
            return tmp


    #Initial parameters
    init_vals = [6e-17,10, redshifts]
    y_fit_guess = model2(wl_fit, *init_vals) + cont


    #Fit
    import scipy.optimize as op
    best_vals, covar = op.curve_fit(model2, wl_fit, flux_fit - cont, sigma=fluxerr_fit, absolute_sigma=True, p0=init_vals)
    print(best_vals)
    z_op = best_vals[-1]
    print("""Curve_fit results:
        Redshift = {0} +- {1} (SDSS: {2})
    """.format(z_op, np.sqrt(covar[-1,-1]), redshifts))

    #Calculate best fit values + confidence values
    y_op = model2(wl_fit, *best_vals) + cont

    from methods import ConfInt
    y_op_lower, y_op_upper = ConfInt(wl_fit, model2, best_vals, covar, [16,84]) + cont

    # Plotting for the paper
    # from plotting import plotting
    # xrng = (10100, 10800)
    # yrng = (3.75e-16, 7.4e-16)
    # fig = plotting.plot_data(wl, flux, xrng=xrng, yrng=yrng, lw=0.8,
    #         label='SDSS0820+1306', label_error='Error',
    #         plot_top_tick=True, z=z_op)
    # ax1, ax2, ax3 = fig.axes
    #
    # ax1.plot(wl_fit, y_op, 'r-', label='[OIII] fit')
    # ax1.fill_between(wl_fit, y_op_lower, y_op_upper, color= 'red', alpha = 0.2)
    #
    # for p in range(len(fit_line_positions)):
    #     xcoord = linelist[p]*(1+z_op)
    #     mask = (wl > xcoord - 1) & (wl < xcoord + 1)
    #     y_val = np.mean(flux[mask])
    #     ax1.axvline(x=xcoord,color='green',linestyle='dashed', lw=1.2)
    #     #ax2.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
    #     ax1.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.15 ),fontsize='small')
    #
    # ax1.legend()
    # import matplotlib
    # x1_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    # ax1.xaxis.set_major_formatter(x1_formatter)
    # y1_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    # ax1.yaxis.set_major_formatter(y1_formatter)
    # for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax3.get_yticklabels() + [ax3.title, ax3.xaxis.label, ax3.yaxis.label] + ax3.get_xticklabels() + ax2.get_yticklabels()):
    #     item.set_fontsize(22)
    # # t = pl.title('SDSS0820')
    # # t.set_y(1.11)
    # #pl.tight_layout()
    #
    # pl.savefig("LineFit.pdf", clobber=True)
    # pl.show(block=True)


    # #Overplot lines
    # for p in range(len(fit_line_positions)):
    #     xcoord = linelist[p]*(1+z_op)
    #     mask = (wl > xcoord - 1) & (wl < xcoord + 1)
    #     y_val = np.mean(flux[mask])
    #     pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
    #     pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.15 ),fontsize='small')
    #
    #
    pl.plot(wl, flux , color = 'black', lw = 0.2, linestyle = 'steps-mid')
    pl.plot(wl, fluxerr, color = 'black', lw = 0.2)
    #pl.plot(wl_fit,flux_fit, color = 'green', lw = 1.0, alpha = 0.5)
    #pl.plot(wl_fit, y_fit_guess)
    pl.plot(wl_fit, y_op, 'r-')
    pl.fill_between(wl_fit, y_op_lower, y_op_upper, color= 'red', alpha = 0.2)
    pl.xlim((10100, 10800))
    pl.ylim((3.75e-16, 7.4e-16))
    pl.xlabel(r'Observed Wavelength  [\AA]')
    pl.ylabel(r'FLux [erg/cm$^2$/s/\AA]')
    pl.tight_layout()
    pl.savefig("Line-fit.pdf", dpi= 150)
    pl.show()
    #Saving to .dat file
    dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64), ("wl sdss", np.float64), ("flux sdss", np.float64) ]
    data = np.array(zip(wl, flux, err_out, wl_sdss, flux_sdss), dtype=dt)
    file_name = "Telluric_corrected_science"
    np.savetxt(i+file_name+".dat", data, header="wl flux fluxerror, sdss_wl, sdss_flux")#, fmt = ['%5.1f', '%2.15E'] )