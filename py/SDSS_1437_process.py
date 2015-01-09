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
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS1437-0147/'
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
    redshifts = 1.30914

    # Load linelist
    fit_line_positions = np.genfromtxt('fitlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)


    from methods import wavelength_conversion
    linelist = wavelength_conversion(linelist, conversion='vacuum_to_air')


    #Cut out fitting region
    mask = np.logical_and(wl > 11350, (wl < 11750))
    wl_fit = wl[mask]
    flux_fit = flux[mask]
    fluxerr_fit = fluxerr[mask]


    #Fit continuum and subtract
    from methods import continuum_fit
    cont, chebfit = continuum_fit(wl_fit, flux_fit, fluxerr_fit, edge_mask_len=20)
    chebfitval = chebyshev.chebval(wl, chebfit)

    pl.plot(wl, flux , color = 'black', lw = 0.2)
    pl.plot(wl, fluxerr, color = 'black', lw = 0.2)
    #pl.show()

    from methods import voigt,gauss

    def model1(t,  amp2, sig22g, sig22l, z):
            tmp = voigt(t, abs(amp2), (1+z)*linelist[3], sig22g, sig22l)
            return tmp

    def model2(t, amp2, sig22g, z):
            tmp = gauss(t, abs(amp2), (1+z)*linelist[3], sig22g)
            return tmp

    init_vals = [6e-17,10, redshifts-0.0002]

    y_fit_guess = model2(wl_fit, *init_vals) + cont

    pl.plot(wl_fit, y_fit_guess)
    # pl.show()

    import scipy.optimize as op

    best_vals, covar = op.curve_fit(model2, wl_fit, flux_fit - cont, sigma=fluxerr_fit, absolute_sigma=True, p0=init_vals)
    print(covar)
    z_op = best_vals[2]
    print("""Curve_fit results:
        Redshift = {0} +- {1} (SDSS: {2})

    """.format(z_op, np.sqrt(covar[2,2]), redshifts))

    y_op = model2(wl, *best_vals)

    pl.plot(wl, y_op + chebfitval, 'r-')

    #Overplot lines
    for p in range(len(fit_line_positions)):
        xcoord = linelist[p]*(1+z_op)
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(flux[mask])
        pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
        pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.4 ),fontsize='x-small')


    pl.plot(wl_fit,flux_fit, color = 'green', lw = 1.0, alpha = 0.5)
    pl.xlim((11000, 11700))
    pl.ylim((0, 2e-15))
    pl.show()