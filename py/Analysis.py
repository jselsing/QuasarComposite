#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2015 Jonatan Selsing"

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')


import numpy as np
from gen_methods import medfilt
import glob

import matplotlib.pylab as pl
import lineid_plot
import seaborn as sns; sns.set_style('ticks')
# cmap = sns.cubehelix_palette(n_colors=6, start=0.0, rot=1.5, gamma=1.0, hue=1.0, light=0.85, dark=0.15, reverse=True, as_cmap=False)
# cmap = sns.color_palette("cubehelix", 3)



def latexify(fig_width=None, fig_height=None, columns=1):
    import matplotlib
    from math import sqrt
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'Qt4Agg',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'font.size': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              # 'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

SPINE_COLOR = 'gray'
def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(True)

    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax


def main():
    latexify()

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

    fig, (ax1, ax2, ax3, ax4) = pl.subplots(ncols=1, nrows=4, sharex=True, figsize=(12, 12))

    #Plotting
    # ax1.plot(wl, power_law2(wl, *popt2) , 'b--')
    # ax1.plot(wl, power_law(wl, *popt) , 'b--')
    ax1.plot(wl, wmean_cont, lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')


    ax2.plot(wl, medfilt(std, 3), lw = 0.5, linestyle = 'steps-mid', label='Normalised variance')

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
    ax1.set_ylim((0.1, 500))


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

    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50])
    ax2.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax2.set_yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2])






    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
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

    for p in range(len(linelist)):
         xcoord = linelist[p]
         mask = (wl > xcoord - 1) & (wl < xcoord + 1)
         y_val = np.mean(wmean_cont[mask])
         # print(y_val)

         ax1.axvline(x=xcoord, ymax=np.log10(y_val)/np.log10(500),  color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax2.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax3.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)
         ax4.axvline(x=xcoord, color='black',linestyle='dashed', lw=0.75, alpha=0.75)

         # trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
         # ann = pl.annotate(fit_line_positions[p,][0], xy=(xcoord, 1.05 ), fontsize='x-small', rotation = 90, xycoords=trans)
         # ax1.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 2.0 ),fontsize='x-small', rotation = 90)


    # fig = pl.gcf()
    # fig.set_size_inches(10,10)
    fig.savefig("../documents/figs/Combined.pdf")
    pl.show()




if __name__ == '__main__':
    main()