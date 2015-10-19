#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function


__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2015 Jonatan Selsing"



from methods import latexify, format_axes, gauss

import numpy as np
from gen_methods import medfilt, smooth
import glob

import matplotlib.pylab as pl
import lineid_plot
import seaborn as sns; sns.set_style('ticks')



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
    std_norm = data_file[:,9]
    wmean_cont = data_file[:,10]




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

    par_guess = [1, -1.70]
    par_guess2 = [1, 5700, -1.7, -1.7]
    wmean[np.where(np.isnan(wmean) == True)] = 0

    mask = (wl > 1300) & (wl < 1350) | (wl > 1425) & (wl < 1475) | (wl > 5500) & (wl < 5800) | (wl > 7300) & (wl < 7500)
    err = ((std)[std != 0])[mask]
    popt, pcov = optimize.curve_fit(power_law, wl[mask], wmean_cont[mask], p0=par_guess, sigma=np.sqrt(err**2 + err_wmean[mask]**2), absolute_sigma=True, maxfev=5000)
    popt2, pcov2 = optimize.curve_fit(power_law2, wl[mask], wmean_cont[mask], p0=par_guess2, sigma=np.sqrt(err**2 + err_wmean[mask]**2), absolute_sigma=True, maxfev=5000)

    print(*popt)
    print(*popt2)

    exit()
    par_guess = [1, -1.7]

    wl_new = wl
    wm = []
    m = []
    geo = []
    med = []
    #Fit
    np.random.seed(12345)
    mask = (wl_new > 1300) & (wl_new < 1350) | (wl_new > 1425) & (wl_new < 1475) | (wl_new > 5500) & (wl_new < 5800) | (wl_new > 7300) & (wl_new < 7500)


    for i in np.arange(10000):
        print('Iteration: ', i)
        err = ((std)[std != 0])[mask]
        resampled_spec = np.random.normal((wmean_cont)[mask], np.sqrt(err**2 + err_wmean[mask]**2))
        popt_wmean, pcov_wmean = optimize.curve_fit(power_law, wl_new[mask], resampled_spec, p0=par_guess,
                                                    sigma=np.sqrt(err**2 + err_wmean[mask]**2), absolute_sigma=True)
        wm.append(popt_wmean)

        err = ((std)[std != 0])[mask]
        resampled_spec = np.random.normal((mean)[mask], np.sqrt(err**2 + err_mean[mask]**2))
        popt_mean, pcov_mean = optimize.curve_fit(power_law, wl_new[mask], resampled_spec, p0=par_guess,
                                                    sigma=np.sqrt(err**2 + err_mean[mask]**2), absolute_sigma=True)
        m.append(popt_mean)

        err = ((std)[std != 0])[mask]
        resampled_spec = np.random.normal((median[std != 0])[mask], err)
        popt_median, pcov_median = optimize.curve_fit(power_law, wl_new[mask], resampled_spec, p0=par_guess,
                                                    sigma=err, absolute_sigma=True, maxfev=600)

        med.append(popt_median)


        err = ((std)[std != 0])[mask]
        resampled_spec = np.random.normal((geo_mean[std != 0])[mask], err)
        # pl.plot(wl_new[mask], resampled_spec)
        popt_geo, pcov_geo = optimize.curve_fit(power_law, wl_new[mask], resampled_spec, p0=par_guess,
                                                    sigma=err, absolute_sigma=True, maxfev=600)


        geo.append(popt_geo)
    # pl.plot(wl, geo_mean)
    # pl.show()


    print("""Composite fit slope wmean...{0} +- {1}""".format(np.mean(wm, axis=0)[1],np.std(wm, axis=0)[1]))
    print("""Composite fit slope mean...{0} +- {1}""".format(np.mean(m, axis=0)[1], np.std(m, axis=0)[1]))
    print("""Composite fit slope median...{0} +- {1}""".format(np.mean(med, axis=0)[1], np.std(med, axis=0)[1]))
    print("""Composite fit slope geo...{0} +- {1}""".format(np.mean(geo, axis=0)[1], np.std(geo, axis=0)[1]))


    # Plotting
    latexify(columns=2, fig_height=7)
    import matplotlib.gridspec as gridspec

    fig = pl.figure()

    gs = gridspec.GridSpec(4, 1, height_ratios=[2,1,1,1])

    ax1 = pl.subplot(gs[0])
    ax2 = pl.subplot(gs[1])
    ax3 = pl.subplot(gs[2])
    ax4 = pl.subplot(gs[3])

    #Plotting
    # ax1.plot(wl, power_law2(wl, *popt2) , 'b--')
    # ax1.plot(wl, power_law(wl, *popt) , 'b--')
    ax1.plot(wl, wmean_cont, lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')

    nbins = len(wl)
    from methods import hist
    log_binned_wl = np.array(hist(wl,[min(wl),max(wl)], int(2*nbins),'log'))
    from scipy.interpolate import InterpolatedUnivariateSpline
    sps = InterpolatedUnivariateSpline(wl, std_norm)
    std_plot = smooth(medfilt(sps(log_binned_wl) , 9), window='hanning', window_len=15)
    wave_std = log_binned_wl

    ax2.plot(wave_std, std_plot, lw = 0.5, linestyle = 'steps-mid', label='Normalised variability')

    ax3.plot(wl, wmean_cont/medfilt(err_wmean, 5), lw = 0.5, linestyle = 'steps-mid', label = 'Signal-to-noise')


    ax4.plot(wl,medfilt( n_spec, 1), label='Number of spectra', lw=0.5)





    #Overplot lines
    fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)

    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(n[1])
        linenames.append(n[0])

    pl.xlabel(r'Rest Wavelength  [$\AA$]')
    ax1.set_ylabel(r'Normalised flux density F$_{\lambda}$')
    ax2.set_ylabel(r'Normalised Variability  $\delta$F$_{\lambda}$')
    ax3.set_ylabel(r'S/N Ratio')
    ax4.set_ylabel(r'Number of spectra')

    ax1.semilogy()
    ax1.semilogx()
    ax1.set_xlim((1000, 11500 ))
    ax1.set_ylim((0.1, 750))


    ax2.semilogy()
    ax2.semilogx()
    ax2.set_xlim((1000, 11500 ))
    ax2.set_ylim((0.001, 90))

    ax3.semilogx()
    ax3.set_xlim((1000, 11500 ))

    ax4.semilogx()
    ax4.set_xlim((1000, 11500 ))
    ax4.set_ylim((0, 9))



    # Formatting axes
    import matplotlib as mpl

    ax4.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax4.set_xticks([1000, 2000, 3000, 5000, 9000])
    ax1.xaxis.set_major_locator(mpl.ticker.NullLocator())
    ax2.xaxis.set_major_locator(mpl.ticker.NullLocator())
    ax3.xaxis.set_major_locator(mpl.ticker.NullLocator())


    ax1.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax2.yaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.set_yticks([0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])
    ax2.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax2.set_yticks([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30])
    ax3.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax3.set_yticks([50, 100, 150, 200, 250, 300])
    ax4.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax4.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])


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
            i.set_linewidth(0.75)

    for p in range(len(linelist)):
         xcoord = linelist[p]
         mask = (wl > xcoord - 1) & (wl < xcoord + 1)
         y_val = np.mean(wmean_cont[mask])
         ax1.vlines(xcoord, ax1.get_ylim()[0], y_val, color='black',linestyle='dashed', lw=0.75, alpha=0.3)
         ax2.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.3)
         ax3.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.3)
         ax4.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.3)

    a = ax1.findobj(mpl.text.Annotation)
    for i in a:
        if '$' in i.get_label():
            i.set_size(10)



    # fig.savefig("../documents/figs/Combined.pdf", rasterized=True, dpi=600)
    pl.show()





if __name__ == '__main__':
    main()