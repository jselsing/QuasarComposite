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

from __future__ import print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"
from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')


import emcee
import triangle
import numpy as np
import matplotlib.pyplot as pl
# import george
# from george import kernels
import glob
from numpy.polynomial import chebyshev



def power_law(params, t):
    a, k = params
    func = a * t ** k
    return func

def power_lawlls(t, a, k):
    func = a * t ** k
    return func


def gauss(params, t):
    amp, mu, sig2 = params
    func = amp * np.exp(-0.5 * (t - mu) ** 2 / sig2)
    return func

def gausslls(t, amp, mu, sig2):
    func = amp * np.exp(-0.5 * (t - mu) ** 2 / sig2)
    return func

def model(params, t):
    func = power_law(params, t)
    return func


def voigtlls(xarr,amp,xcen,Gfwhm,Lfwhm):
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



def lnprior_base(p):
    a, k = p
    if not 0 < a < 20:
        return -np.inf
    if not -5 < k < 5:
        return -np.inf

    return 0.0

def lnlike_gp(p, t, y, yerr):
    a, tau = np.exp(p[:2])
    gp = george.GP(a * kernels.Matern32Kernel(tau))
    gp.compute(t, yerr)
    return gp.lnlikelihood(y - model(p[2:], t))


def lnprior_gp(p):
    lna, lntau = p[:2]
    if not -5 < lna < 5:
        return -np.inf
    if not -5 < lntau < 5:
        return -np.inf
    return lnprior_base(p[2:])


def lnprob_gp(p, t, y, yerr):
    lp = lnprior_gp(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_gp(p, t, y, yerr)


def fit_gp(initial, data, nwalkers=32):
    ndim = len(initial)
    p0 = [np.array(initial) + 1e-1 * np.random.randn(ndim)
          for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, args=data)

    print("Running burn-in")
    p0, lnp, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print("Running second burn-in")
    p = p0[np.argmax(lnp)]
    p0 = [p + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    p0, _, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print("Running production")
    p0, _, _ = sampler.run_mcmc(p0, 1000)

    return sampler

def generate_data(params, N, rng=(0, 15)):
    gp = george.GP(0.00 * kernels.ExpSquaredKernel(0.01))
    t = rng[0] + np.diff(rng) * np.sort(np.random.rand(N))
    y = gp.sample(t)
    print(params)
    y += model(params, t)
    yerr = 0.05 + 0.05 * np.random.rand(N)
    y += yerr * np.random.randn(N)
    return t, y, yerr



def test():
    """ Testing Docstring"""
    pass


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
    return chebfitval


if __name__ == '__main__':
    #Files

    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/')
    redshifts = [1.1257, 1.98041, 1.57697, 1.82389, 1.51237, 2.096, 1.30914]
    print(np.mean(redshifts))
    for x,i in enumerate(sdssobjects):

        obs = np.genfromtxt(glob.glob(i+'Telluric_corrected_science.dat')[0])


        fit_line_positions = np.genfromtxt('fitlinelist.txt', dtype=None)

        wl = obs[:,0]
        flux = obs[:,1]
        fluxerr = obs[:,2]
        sdss_wl = (obs[:,3])[np.where(obs[:,3] != 0)]
        sdss_flux = (obs[:,4])[np.where(obs[:,3] != 0)]


        linelist = []
        for n in fit_line_positions:
            linelist.append(n[1])
        linelist = np.array(linelist)
        linelist = linelist / (1.0 + 2.735182E-4 + 131.4182 / linelist**2 + 2.76249E8 / linelist**4)

        fac = 2*redshifts[x]
        reg = 200 * fac
        inner = 4 * fac
        # xcoord1 = linelist[0]*(1+redshifts[x])
        # mask1 = (wl > xcoord1 - reg) & (wl < xcoord1 + reg)

        xcoord2 = linelist[0]*(1+redshifts[x])
        mask2 = (wl > xcoord2 - reg) & (wl < xcoord2 + reg)

        xcoord3 = linelist[1]*(1+redshifts[x])
        mask3 = (wl > xcoord3 - reg) & (wl < xcoord3 + reg)

        xcoord4 = linelist[2]*(1+redshifts[x])
        mask4 = (wl > xcoord4 - reg) & (wl < xcoord4 + reg)


        mask = mask2 & mask3 & mask4
        # mask = mask1 & mask2

        wl_fit = wl[mask]
        flux_fit = flux[mask]
        fluxerr_fit = fluxerr[mask]



        # wlrng = max(wl_fit) - min(wl_fit)
        # mask_cont = np.logical_and(np.logical_or((wl < np.mean(wl_fit) - wlrng/4),(wl > np.mean(wl_fit) + wlrng/4)), (wl < 21000))
        # chebfit = chebyshev.chebfit(wl[mask_cont], flux[mask_cont], deg = 5, w=1/fluxerr[mask_cont]**2)
        # chebfitval = chebyshev.chebval(wl, chebfit)
        # pl.plot(wl, chebfitval)
        # chebfitval = chebyshev.chebval(wl_fit, chebfit)
        chebfitval = continuum_fit(wl_fit, flux_fit, fluxerr_fit)
        #=============================================================----------------------------------------------------------
        # Least squares
        #=============================================================----------------------------------------------------------
        # Do the least-squares fit and compute the uncertainties.



        # ## Mg II
        # def modelllsmgii(t, amp, mu, amp2, mu2, sig2):
        #     func = gausslls(t, amp, mu, sig2) + gausslls(t, amp2, mu2, sig2)
        #     return func
        # import scipy.optimize as op
        # init_vals = [1.5e-16, linelist[0]*(1+redshifts[x]), 1.5e-16, linelist*(1+redshifts[x]), 50]
        # best_vals, covar = op.curve_fit(modelllsmgii, wl_fit, flux_fit - chebfitval, sigma=fluxerr_fit - chebfitval, absolute_sigma=True, p0=init_vals)
        # amp_op, mu_op, amp2_op, mu2_op, sig2_op = best_vals
        # print("""Curve_fit results:
        #     amp = {0} +- {1} (Guess: {2})
        #     mu = {3} +- {4} (Guess: {5})
        #     amp2 = {6} +- {7} (Guess: {8})
        #     mu2 = {9} +- {10} (Guess: {11})
        #     sigma**2 = {12} +- {13} (Guess: {14})
        # """.format(amp_op, np.sqrt(covar[0, 0]), init_vals[0], mu_op, np.sqrt(covar[1, 1]), init_vals[1],
        #            amp2_op, np.sqrt(covar[2, 2]), init_vals[2], mu2_op, np.sqrt(covar[3, 3]), init_vals[3],
        #            sig2_op, np.sqrt(covar[4, 4]), init_vals[4]))
        # # Plot the least-squares result.
        # y_op_tot = modelllsmgii(wl_fit, amp_op, mu_op, amp2_op, mu2_op, sig2_op) + chebfitval
        # y_op_comp1 = modelllsmgii(wl_fit, amp_op, mu_op, 0, mu2_op, sig2_op) + chebfitval
        # y_op_comp2 = modelllsmgii(wl_fit, 0, mu_op, amp2_op, mu2_op, sig2_op) + chebfitval
        # pl.plot(wl_fit, y_op_tot, "r")
        # pl.plot(wl_fit, y_op_comp1, "--b")
        # pl.plot(wl_fit, y_op_comp2, "--b")

        ## Hbeta + [OIII]
        def modelllshboiii(t, amp, mu, sig2g, sig2l, amp2, sig22g, sig22l, amp3, amp4, z):
            func = voigtlls(t, abs(amp), mu, sig2g, sig2l) + voigtlls(t, abs(amp2), (1+z)*linelist[1], sig22g, sig22l) + voigtlls(t, abs(amp3),
                            (1+z)*linelist[2], sig22g, sig22l) + voigtlls(t, abs(amp4), (1+z)*linelist[3], sig22g, sig22l)
            return func



        import scipy.optimize as op
        init_vals = [2.4e-16, linelist[1]*(1+redshifts[x]), 60, 1,
                     0.3e-16,
                     6, 1,
                     0.3e-16,
                     1e-16, redshifts[x]]
        from gen_methods import smooth
#smooth(flux_fit, window_len=10, window='hanning')
        best_vals, covar = op.curve_fit(modelllshboiii, wl_fit, flux_fit - chebfitval, sigma=fluxerr_fit, absolute_sigma=True, p0=init_vals, maxfev=100000)
        amp_op, mu_op, sig2g_op, sig2l_op, amp2_op, sig22g_op, sig22l_op, amp3_op, amp4_op, z_op = best_vals

        print("""Curve_fit results:
            Redshift = {0} +- {1} (SDSS: {2})

        """.format(z_op, np.sqrt(covar[9,9]), redshifts[x]))
        print(amp_op, mu_op, sig2g_op, sig2l_op, amp2_op, sig22g_op, sig22l_op, amp3_op, amp4_op, z_op)
        # Plot the least-squares result.
        y_op_tot= modelllshboiii(wl_fit, amp_op, mu_op, sig2g_op, sig2l_op, amp2_op, sig22g_op, sig22l_op, amp3_op, amp4_op, z_op) + chebfitval
        y_op_comp1 = modelllshboiii(wl_fit,amp_op, mu_op, sig2g_op, sig2l_op, 0, sig22g_op, sig22l_op, 0, 0, z_op) + chebfitval
        y_op_comp2 = modelllshboiii(wl_fit, 0, mu_op, sig2g_op, sig2l_op, amp2_op, sig22g_op, sig22l_op, 0, 0, z_op) + chebfitval
        y_op_comp3 = modelllshboiii(wl_fit, 0, mu_op, sig2g_op, sig2l_op, 0, sig22g_op, sig22l_op, amp3_op, 0, z_op) + chebfitval
        y_op_comp4 = modelllshboiii(wl_fit, 0, mu_op, sig2g_op, sig2l_op, 0, sig22g_op, sig22l_op, 0, amp4_op, z_op) + chebfitval
        pl.plot(wl_fit, y_op_tot, "r")
        pl.plot(wl_fit, y_op_comp1, "--b")
        pl.plot(wl_fit, y_op_comp2, "--b")
        pl.plot(wl_fit, y_op_comp3, "--b")
        pl.plot(wl_fit, y_op_comp4, "--b")









        from gen_methods import smooth
        #smooth(flux_fit, window_len=10, window='hanning')
        # pl.plot(wl,smooth(flux, window_len=1, window='hanning'), lw=0.5, color = 'black', label='Corrected')
        from xshoo.binning import binning1d
        pl.plot(binning1d(wl, 5), smooth(binning1d(flux, 5, err= fluxerr)[0], window_len=2, window='hanning') , color = 'black', lw = 0.5, alpha=0.7)
        pl.xlabel(r'Observed Wavelength  [\AA]')
        pl.ylabel(r'FLux Error [erg/cm$^2$/s/\AA]')
        pl.title('Line-fit')
        #=============================================================----------------------------------------------------------

        #=============================================================----------------------------------------------------------
        # Maximum Likelihood
        #=============================================================----------------------------------------------------------
        # # Find the maximum likelihood value.
        # chi2 = lambda *args: -2 * lnlike_gp(*args)
        # init_vals = [0.0, 0.0] + [1e-9, -1.6]
        # result = op.minimize(chi2, [init_vals], args=(wl_out, flux_out, err_out))
        #
        # a_ml, b_ml, a2_ml, b2_ml = result["x"]
        # covar_ml = result["hess_inv"]
        # print("""Maximum likelihood result:
        #     a = {0} (Guess: {1})
        #     b = {2} (Guess: {3})
        #     a2 = {4} (Guess: {5})
        #     b2 = {6} (Guess: {7})
        #
        # """.format(a_ml, init_vals[0], b_ml, init_vals[1],
        #            a2_ml, init_vals[2], b2_ml, init_vals[3]))
        #
        # # Plot the least-squares result.
        # y_ml = modellls(wl_out, a_ml, b_ml, a2_ml, b2_ml)
        # pl.plot(wl_out, y_ml, "--g")
        pl.xlim((0.99*min(wl[mask]), 1.01*max(wl[mask])))
        pl.ylim((0.7*np.median(flux[mask]), 2.0* np.median(flux[mask])))

        composite = np.genfromtxt('linelist.txt', dtype=None)
        linelist = []
        for n in composite:
            linelist.append(n[1])
        linelist = np.array(linelist)
        linelist = linelist / (1.0 + 2.735182E-4 + 131.4182 / linelist**2 + 2.76249E8 / linelist**4)
        for p in range(len(composite)):
            xcoord = linelist[p]*(1+redshifts[x])
            mask = (wl > xcoord - 1) & (wl < xcoord + 1)
            y_val = np.median(flux[mask])
            pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.5)
            pl.annotate(composite[p,][0],xy=(xcoord, y_val * 1.2 ),fontsize='x-small')

        pl.show(block=True)






