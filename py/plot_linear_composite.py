#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


from methods import latexify, format_axes, gauss

import numpy as np
from gen_methods import medfilt
import glob

import matplotlib.pylab as pl
import lineid_plot
import seaborn as sns; sns.set_style('ticks')
# cmap = sns.cubehelix_palette(n_colors=6, start=0.0, rot=1.5, gamma=1.0, hue=1.0, light=0.85, dark=0.15, reverse=True, as_cmap=False)
cmap = sns.color_palette("cubehelix", 3)




def main():
    # latexify()
    import numpy as np

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

    from scipy.interpolate import InterpolatedUnivariateSpline
    mask = (np.where(err_wmean != 0))
    f = InterpolatedUnivariateSpline(wl[mask], wmean_cont[mask], w=err_wmean[mask], k=5)
    wmean_cont = f(wl)
    f = InterpolatedUnivariateSpline(wl[mask], err_wmean[mask], k=5)
    err_wmean = f(wl)




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


    par_guess = [1, -1.7]
    par_guess2 = [1, 5700, -1.7, -1.7]
    wmean[np.where(np.isnan(wmean) == True)] = 0
    mask = (wl > 1300) & (wl < 1350) | (wl > 1425) & (wl < 1475) | (wl > 5500) & (wl < 5800) | (wl > 7300) & (wl < 7500) #| (wl > 9700) & (wl < 9900) | (wl > 10200) & (wl < 10600)
    popt, pcov = optimize.curve_fit(power_law, wl[mask], wmean_cont[mask], p0=par_guess, sigma=err_wmean[mask], absolute_sigma=True)
    popt2, pcov2 = optimize.curve_fit(power_law2, wl[mask], wmean_cont[mask], p0=par_guess2, sigma=err_wmean[mask], absolute_sigma=True)

    print(*popt)
    print(*popt2)



    #Plotting
    latexify(columns=2)
    fig, ax = pl.subplots()

    ax.plot(wl, medfilt(wmean_cont, 5),
            lw = 0.5, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite', color=cmap[1])

    ax.plot(wl, power_law(wl, *popt),
            linestyle='dashed', label ='Pure power law fit', color=cmap[2])

    sdss_compo = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/sdss_compo.dat')
    sdss_wl = sdss_compo[:,0]
    sdss_flux = sdss_compo[:, 1]

    norm_reg = 1430

    mask = (wl > norm_reg) & (wl < norm_reg + 20)
    norm1 = np.median(wmean_cont[mask])
    mask = (sdss_wl > norm_reg) & (sdss_wl < norm_reg + 20)
    norm2 = np.median(sdss_flux[mask])


    ax.plot(sdss_wl, sdss_flux * (norm1/norm2),
            linestyle='solid', label ='Full sample SDSS composite', color=cmap[0])


    #Overplot lines
    fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)

    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(n[1])
        linenames.append(n[0])

    pl.semilogy()

    # Formatting axes
    import matplotlib as mpl

    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.get_xaxis().tick_bottom()
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_yticks([0.3, 1, 3, 10, 30, 100, 300])
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())

    # pl.legend(loc=3)
    ax.set_xlim((1000, 11500))
    ax.set_ylim((0.2, 500))
    format_axes(ax)





    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel(r'Rescaled flux density F$_\lambda$')
    pl.tight_layout()




    val = []
    for p in range(len(linelist)):
        xcoord = linelist[p]
        mask = (wl > xcoord - 1) & (wl < xcoord + 1)
        y_val = np.mean(wmean_cont[mask])
        val.append(2 * y_val)
    print(val)
    arrow_tips = val
    lineid_plot.plot_line_ids(wl, wmean_cont, linelist, linenames, arrow_tip=arrow_tips, ax=ax)
    for i in ax.lines:
        if '$' in i.get_label():
            i.set_alpha(0.3)
    a = ax.findobj(mpl.text.Annotation)
    for i in a:
        if '$' in i.get_label():
            i.set_size(10)


    pl.savefig('../documents/figs/compo_full_sample.pdf')
    pl.show()



if __name__ == '__main__':
    main()