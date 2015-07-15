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
# cmap = sns.cubehelix_palette(n_colors=6, start=0.0, rot=1.5, gamma=1.0, hue=1.0, light=0.85, dark=0.15, reverse=True, as_cmap=False)
# cmap = sns.color_palette("cubehelix", 3)





def main():


    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'


    data_file = np.genfromtxt(root_dir+'Composite.dat')

    wl = data_file[:,0]
    mean = data_file[:,1]
    err_mean = data_file[:,2]
    wmean = data_file[:,3]
    err_wmean = data_file[:,4]
    geo_mean = data_file[:,5]
    median = data_file[:,6]
    n_spec = data_file[:,7]
    std = data_file[:,8]
    # wl_sdss = data_file[:,9]
    # mean_sdss = data_file[:,10]
    wmean_cont = data_file[:,9]



    #Fitting power laws
    from scipy import optimize

    def power_law(x_tmp, a_tmp, k_tmp):

        tmp = a_tmp * x_tmp ** k_tmp
        return tmp

    def power_law2(x_tmp, a1_tmp, x_c, k1_tmp, k2_tmp):

        tmp1 = power_law(x_tmp, a1_tmp,k1_tmp)[x_tmp<x_c]
        scale2loc = np.argmin(np.abs(x_tmp - x_c))
        a2_tmp = power_law(x_tmp[scale2loc], a1_tmp, (k1_tmp - k2_tmp))

        tmp2 = power_law(x_tmp, a2_tmp,k2_tmp)[x_tmp>= x_c]

        return np.concatenate((tmp1,tmp2))

    def power_law3(x_tmp, a_tmp, k_tmp, b_tmp):

        tmp = a_tmp * x_tmp ** (k_tmp + b_tmp * x_tmp)
        return tmp

    par_guess = [1, -1.7]
    par_guess2 = [1, 5700, -1.7, -1.7]
    wmean[np.where(np.isnan(wmean) == True)] = 0

    mask = (wl > 1300) & (wl < 1350) | (wl > 1425) & (wl < 1475) | (wl > 5500) & (wl < 5800) | (wl > 7300) & (wl < 7500)

    popt, pcov = optimize.curve_fit(power_law, wl[mask], wmean_cont[mask], p0=par_guess, sigma=err_wmean[mask], absolute_sigma=True)
    popt2, pcov2 = optimize.curve_fit(power_law2, wl[mask], wmean_cont[mask], p0=par_guess2, sigma=err_wmean[mask], absolute_sigma=True)

    print(*popt)
    print(*popt2)

    #Plotting
    ratio = (1.0 + np.sqrt(5.0))/2.0
    y_size = 12
    x_size = 12

    import matplotlib.gridspec as gridspec

    fig = pl.figure(figsize=(x_size, y_size))

    gs = gridspec.GridSpec(4, 1, height_ratios=[2,1,1,1])

    ax1 = pl.subplot(gs[0])
    ax2 = pl.subplot(gs[1])
    ax3 = pl.subplot(gs[2])
    ax4 = pl.subplot(gs[3])

    # fig, (ax1, ax2, ax3, ax4) = pl.subplots(ncols=1, nrows=4, sharex=True, figsize=(12, 12))

    #Plotting
    # ax1.plot(wl, power_law2(wl, *popt2) , 'b--')
    # ax1.plot(wl, power_law(wl, *popt) , 'b--')
    ax1.plot(wl, wmean_cont, lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')

    nbins = len(wl)
    from methods import hist
    log_binned_wl = np.array(hist(wl,[min(wl),max(wl)], int(2*nbins),'log'))
    from scipy.interpolate import InterpolatedUnivariateSpline
    sps = InterpolatedUnivariateSpline(wl, std)
    std_plot = smooth(medfilt(sps(log_binned_wl) , 9), window='hanning', window_len=15)
    wave_std = log_binned_wl

    ax2.plot(wave_std, std_plot, lw = 0.5, linestyle = 'steps-mid', label='Normalised variance')

    ax3.plot(wl, wmean_cont/medfilt(err_wmean, 5), lw = 0.5, linestyle = 'steps-mid', label = 'Signal-to-noise')


    ax4.plot(wl,medfilt( n_spec, 151), label='Number of spectra')





    #Overplot lines
    fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)

    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(n[1])
        linenames.append(n[0])





    pl.xlabel(r'Rest Wavelength  [$\AA$]')
    ax1.set_ylabel(r'Normalised flux density F_{\lambda}')
    ax2.set_ylabel(r'Normalised Variance')

    ax3.set_ylabel(r'S/N Ratio')
    ax4.set_ylabel(r'Number of spectra')




    ax1.semilogy()
    ax1.semilogx()
    ax1.set_xlim((1000, 11500 ))
    ax1.set_ylim((0.05, 500))


    ax2.semilogy()
    ax2.semilogx()
    ax2.set_xlim((1000, 11500 ))
    ax2.set_ylim((1e-2, 2.0))

    ax3.semilogx()
    ax3.set_xlim((1000, 11500 ))

    ax4.semilogx()
    ax4.set_xlim((1000, 11500 ))
    ax4.set_ylim((0, 9))



    # Formatting axes
    import matplotlib as mpl

    ax4.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax4.set_xticks([1000, 2000, 3000, 5000, 9000])
    # ax4.get_xaxis().tick_bottom()
    # ax4.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax1.xaxis.set_major_locator(mpl.ticker.NullLocator())
    ax2.xaxis.set_major_locator(mpl.ticker.NullLocator())
    ax3.xaxis.set_major_locator(mpl.ticker.NullLocator())


    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200])
    ax2.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax2.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5, 1])
    ax3.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax3.set_yticks([20, 40, 60, 80, 100])
    ax4.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax4.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])




    # ax1.legend()
    # ax2.legend()
    # ax3.legend()
    # ax4.legend()
    pl.tight_layout()
    fig.subplots_adjust(hspace=0)

    format_axes(ax1)
    format_axes(ax2)
    format_axes(ax3)
    format_axes(ax4)



    val = []
    for p in range(len(linelist)):
        xcoord = linelist[p]
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(wmean_cont[mask])
        val.append(1.5 * y_val)
    arrow_tips = val
    lineid_plot.plot_line_ids(wl, wmean_cont, linelist, linenames, arrow_tip=arrow_tips, ax=ax1)
    for i in ax1.lines:
        if '$' in i.get_label():
            i.set_alpha(0.3)

    xl = ax1.get_ylim()[1] - ax1.get_ylim()[0]

    for p in range(len(linelist)):
         xcoord = linelist[p]
         mask = (wl > xcoord - 1) & (wl < xcoord + 1)
         y_val = np.mean(wmean_cont[mask])
         print((y_val/np.log10(xl)))

         ax1.vlines(xcoord, ax1.get_ylim()[0], y_val, color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax2.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax3.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax4.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)

         # trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
         # ann = pl.annotate(fit_line_positions[p,][0], xy=(xcoord, 1.05 ), fontsize='x-small', rotation = 90, xycoords=trans)
         # ax1.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 2.0 ),fontsize='x-small', rotation = 90)


    # fig = pl.gcf()
    # fig.set_size_inches(10,10)
    latexify(fig_width=x_size, fig_height=y_size)
    fig.savefig("../documents/figs/Combined.pdf")
    pl.show()




if __name__ == '__main__':
    main()