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

def test():
    """ Testing Docstring"""
    pass

def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
    from scipy.special import wofz

    """
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Converted from
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    """
    tmp = 1.0/wofz(np.zeros((len(xarr)))+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    tmp = tmp*amp*wofz(2*np.sqrt(np.log(2.0))*(xarr-xcen)/Gfwhm+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    return tmp


def continuum_fit(wl, flux, fluxerr, edge_mask_len = 300):
    """
    Small function to estimate continuum. Takes spectrum and uses only the edges, as specified by edge_mask_len, fits
    a polynomial and returns the fit.
    :param wl: Wavelenght array in which to estimate continuum
    :param flux: Corresponding flux array
    :param fluxerr: Corresponding array containing flux errors
    :param edge_mask_len: Size of edges to use to interpolate continuum
    :return: Interpolated continuum
    """
    wl_chebfit = np.hstack((wl[:edge_mask_len],wl[-edge_mask_len:]))
    mask_cont = np.array([i for i,n in enumerate(wl) if n in wl_chebfit])
    chebfit = chebyshev.chebfit(wl_chebfit, flux[mask_cont], deg = 1, w=1/fluxerr[mask_cont]**2)
    chebfitval = chebyshev.chebval(wl, chebfit)
    return chebfitval, chebfit








if __name__ == '__main__':

    #Get spectrum
    object = []
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS1150-0023/'
    obs = np.genfromtxt(glob.glob(root_dir+'Telluric_corrected_science.dat')[0])
    wl = obs[:,0]
    flux = obs[:,1]
    fluxerr = obs[:,2]
    sdss_wl = (obs[:,3])[np.where(obs[:,3] != 0)]
    sdss_flux = (obs[:,4])[np.where(obs[:,3] != 0)]
    redshifts = 1.98041


    # Load linelist
    fit_line_positions = np.genfromtxt('fitlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)
    linelist = linelist / (1.0 + 2.735182E-4 + 131.4182 / linelist**2 + 2.76249E8 / linelist**4)



    #Cut out fitting region
    mask = np.logical_and(wl > 14800, (wl < 15000))
    wl_fit = wl[mask]
    flux_fit = flux[mask]
    fluxerr_fit = fluxerr[mask]

    #Fit continuum and subtract
    cont, chebfit = continuum_fit(wl_fit, flux_fit, fluxerr_fit, edge_mask_len=2)
    chebfitval = chebyshev.chebval(wl, chebfit)

    #Define models to use
    from methods import voigt,gauss
    def model1(t,  amp2, sig22g, sig22l, z):
            tmp = voigt(t, abs(amp2), (1+z)*linelist[2], sig22g, sig22l)
            return tmp

    def model2(t, amp2, sig22g, z):
            tmp = gauss(t, abs(amp2), (1+z)*linelist[2], sig22g)
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
    y_op = model2(wl, *best_vals) + chebfitval

    from methods import ConfInt
    y_op_lower, y_op_upper = ConfInt(wl_fit, model2, best_vals, covar, [16,84]) + cont

    #Overplot lines
    for p in range(len(fit_line_positions)):
        xcoord = linelist[p]*(1+z_op)
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(flux[mask])
        pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
        pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.4 ),fontsize='x-small')


    pl.plot(wl, flux , color = 'black', lw = 0.2, linestyle = 'steps-mid')
    pl.plot(wl, fluxerr, color = 'black', lw = 0.2)
    #pl.plot(wl_fit,flux_fit, color = 'green', lw = 1.0, alpha = 0.5)
    #pl.plot(wl_fit, y_fit_guess)
    pl.plot(wl, y_op, 'r-')
    pl.fill_between(wl_fit, y_op_lower, y_op_upper, color= 'red', alpha = 0.2)
    pl.xlim((14000, 15300))
    pl.ylim((1e-17, 5e-16))
    pl.show()