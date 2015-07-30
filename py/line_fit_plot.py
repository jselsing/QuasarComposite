#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2015 Jonatan Selsing"


from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')

from methods import latexify, format_axes, gauss

import numpy as np
from gen_methods import medfilt, smooth
import glob

import matplotlib.pylab as pl
import lineid_plot
import seaborn as sns; sns.set_style('ticks')
cmap = sns.color_palette("cubehelix", 4)
# def model2(t, amp2, sig22g, z):
#         tmp = gauss(t, abs(amp2), (1+z)*linelist[2], sig22g)
#         return tmp



if __name__ == '__main__':



    # datfile = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS0820+1306/Telluric_corrected_science.dat')
    # datfile = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS1437-0147/Telluric_corrected_science.dat')
    datfile = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS1354-0013/Telluric_corrected_science.dat')

    # z = 1.124
    z = 1.5124

    norm = 1e16

    # wl = datfile[:,0]
    # bp_map = datfile[:,3]
    # flux = np.ma.array(datfile[:,1], mask=bp_map)
    # fluxerror = np.ma.array(datfile[:,2], mask=bp_map)

    wl = datfile[:,0]/ (1 + z)
    bp_map = datfile[:,3]
    flux = np.ma.array(datfile[:,1]*(1 + z)*norm, mask=bp_map)
    fluxerror = np.ma.array(datfile[:,2]*(1 + z)*norm, mask=bp_map)

 # Load linelist
    fit_line_positions = np.genfromtxt('data/fitlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)


    from methods import wavelength_conversion
    linelist = wavelength_conversion(linelist, conversion='vacuum_to_air')

    #Cut out fitting region
    mask = ((wl >= 4975) & (wl <= 5100)) & (flux.mask == 0.)
    wl_fit = wl[mask]
    flux_fit = flux[mask]
    fluxerr_fit = fluxerror[mask]

    fluxerr_new = []
    for j, (k, l) in enumerate(zip(flux_fit,fluxerr_fit)):
        if k > 1.5 * flux_fit[j-2] and k > 0:
            fluxerr_new.append(l*50)
        elif k < 0.75 * flux_fit[j-2] and k > 0:
            fluxerr_new.append(l*50)
        else:
            fluxerr_new.append(l)
    from gen_methods import smooth
    fluxerr_fit = smooth(np.array(fluxerr_new), window_len=15, window='hanning')


    #Fit continuum and subtract
    from methods import continuum_fit
    from numpy.polynomial import chebyshev
    cont, chebfit = continuum_fit(wl_fit, flux_fit, fluxerr_fit, edge_mask_len=20)
    chebfitval = chebyshev.chebval(wl, chebfit)


    #Define models to use
    from methods import gauss

    def model2(t, amp2, sig22g, z):
            tmp = gauss(t, abs(amp2), (1+z)*linelist[2], sig22g)
            return tmp


    #Initial parameters
    init_vals = [1e-15*norm,10, 0]
    y_fit_guess = model2(wl_fit, *init_vals) + cont

    #Fit
    import scipy.optimize as op
    np.random.seed(12345)
    y_op = []
    vals = []
    for i in np.arange(130):
        # print('Iteration: ', i)
        resampled_spec = np.random.normal(flux_fit, fluxerr_fit)

        cont, chebfit = continuum_fit(wl_fit, resampled_spec, fluxerr_fit, edge_mask_len=20)
        chebfitval = chebyshev.chebval(wl, chebfit)

        best_vals, covar = op.curve_fit(model2, wl_fit, resampled_spec - cont, sigma=fluxerr_fit, absolute_sigma=True, p0=init_vals)
        # print(best_vals)




        #Calculate best fit values + confidence values
        y_op.append(model2(wl_fit, *best_vals) + cont)
        vals.append(best_vals)


    print(np.mean(vals, axis = 0))

    y_op_lower, y_op_upper = np.mean(y_op, axis=0) - 2*np.std(y_op, axis=0), np.mean(y_op, axis=0) + 2*np.std(y_op, axis=0)
    y_op = np.mean(y_op, axis=0)
    # def autocorr(x):
    #     result = np.correlate(x, x, mode='full')
    #     return result[result.size/2:]
    #
    # pl.plot(autocorr(flux)/np.median(autocorr(flux)))
    # pl.plot((flux/np.median((flux))))
    # pl.show()
    # from methods import ConfInt
    # y_op_lower, y_op_upper = ConfInt(wl_fit, model2, best_vals, covar, [16,84]) + cont


    #Plotting
    ratio = (1.0 + np.sqrt(5.0))/2.0
    latexify()
    fig, ax = pl.subplots()
    # latexify(fig_width=4*ratio, fig_height=4, columns=1)
    # latexify(columns=1)
    # fig, ax = pl.subplots()

    from gen_methods import medfilt
    # ax.plot(wl, medfilt(flux, 3)/norm, drawstyle='steps-mid', lw = 0.5, alpha=0.5)
    ax.errorbar(wl, medfilt(flux, 3)/norm, yerr=fluxerror/norm, fmt='o', marker='o', capsize=0.75, ms=0.25, mew=0.25, color='black', elinewidth=0.15, alpha=0.2)


    ax.plot(wl_fit, y_op/norm, lw=1, color = cmap[2])
    ax.fill_between(wl_fit, y_op_lower/norm, y_op_upper/norm, alpha = 0.5, color = cmap[2])




    ax.set_xlim((4700, 5100))
    ax.set_ylim((5e-16, 9e-16))




    # ax.set_xticks([4750, 4850, 4950, 5050, 5150])
    format_axes(ax)
    ax.set_xlabel(r'Restframe Wavelength [$\AA$]')
    ax.set_ylabel(r'Flux density [erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$]')
    pl.tight_layout()
    import matplotlib as mpl
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xticks([4700, 4800, 4900, 5000, 5100])
    #Overplot lines
    fit_line_positions = np.genfromtxt('data/lineplotlinelist.txt', dtype=None)
    import lineid_plot

    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(n[1])
        linenames.append(n[0])
    linelist = np.array(linelist)
    linelist = wavelength_conversion(linelist, conversion='vacuum_to_air')


    lineid_plot.plot_line_ids(wl, medfilt(flux, 3)/norm, linelist, linenames, ax=ax)
    for i in ax.lines:
        if '$' in i.get_label():
            i.set_alpha(0.3)
            i.set_linewidth(0.75)




    fig.tight_layout()





    pl.savefig("../documents/figs/LineFit.pdf", clobber=True)
    pl.show()