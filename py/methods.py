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

__all__ = ["voigt", "gauss", "continuum_fit", "wavelength_conversion"]
__version__ = "0.0.1"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2015 Jonatan Selsing"

import numpy as np


def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
    """
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Converted from
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    """
    from scipy.special import wofz
    tmp = 1.0/wofz(np.zeros((len(xarr)))+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    tmp = tmp*amp*wofz(2*np.sqrt(np.log(2.0))*(xarr-xcen)/Gfwhm+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    return tmp

def gauss(xarr, amp, mu, sig2):
    """
    Generic gaussian function
    :param xarr: x-array
    :param amp: amplitude
    :param mu: x-center
    :param sig2: variance
    :return: gaussian function
    """
    tmp = amp * np.exp(-0.5 * (xarr - mu) ** 2 / sig2)
    return tmp


def continuum_fit(wl, flux, fluxerr, edge_mask_len=50):
    """
    Small function to estimate continuum. Takes spectrum and uses only the edges, as specified by edge_mask_len, fits
    a polynomial and returns the fit.
    :param wl: Wavelenght array in which to estimate continuum
    :param flux: Corresponding flux array
    :param fluxerr: Corresponding array containing flux errors
    :param edge_mask_len: Size of edges to use to interpolate continuum
    :return: Interpolated continuum
    """
    from numpy.polynomial import chebyshev
    wl_chebfit = np.hstack((wl[:edge_mask_len],wl[-edge_mask_len:]))
    mask_cont = np.array([i for i,n in enumerate(wl) if n in wl_chebfit])
    chebfit = chebyshev.chebfit(wl_chebfit, flux[mask_cont], deg = 1, w=1/fluxerr[mask_cont]**2)
    chebfitval = chebyshev.chebval(wl, chebfit)
    return chebfitval, chebfit

def wavelength_conversion(input_wavelenght, conversion = 'vacuum_to_air'):
    """Convert wavelenghts between air and vacuum.
    Coversion from vacuum to air, taken from: https://www.sdss3.org/dr9/spectro/spectro_basics.php"""
    if conversion == 'vacuum_to_air':
        output_wavelenght = input_wavelenght / (1.0 + 2.735182E-4 + 131.4182 / input_wavelenght**2 + 2.76249E8 / input_wavelenght**4)

    return output_wavelenght


def test():
    """ Testing Docstring"""
    pass


if __name__ == '__main__':
    test()