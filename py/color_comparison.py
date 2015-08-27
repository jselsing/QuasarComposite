#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')

from methods import latexify, format_axes, gauss

def _set_spine_position(spine, position):
    """
    Set the spine's position without resetting an associated axis.

    As of matplotlib v. 1.0.0, if a spine has an associated axis, then
    spine.set_position() calls axis.cla(), which resets locators, formatters,
    etc.  We temporarily replace that call with axis.reset_ticks(), which is
    sufficient for our purposes.
    """
    axis = spine.axis
    if axis is not None:
        cla = axis.cla
        axis.cla = axis.reset_ticks
    spine.set_position(position)
    if axis is not None:
        axis.cla = cla

def despine(fig=None, ax=None, top=True, right=True, left=False,
            bottom=False, offset=None, trim=False):
    """Remove the top and right spines from plot(s).

    fig : matplotlib figure, optional
        Figure to despine all axes of, default uses current figure.
    ax : matplotlib axes, optional
        Specific axes object to despine.
    top, right, left, bottom : boolean, optional
        If True, remove that spine.
    offset : int or None  (default), optional
        Absolute distance, in points, spines should be moved away
        from the axes (negative values move spines inward).
    trim : bool, optional
        If true, limit spines to the smallest and largest major tick
        on each non-despined axis.

    Returns
    -------
    None

    """
    # Get references to the axes we want
    if fig is None and ax is None:
        axes = pl.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        for side in ["top", "right", "left", "bottom"]:
            # Toggle the spine objects
            is_visible = not locals()[side]
            ax_i.spines[side].set_visible(is_visible)
            if offset is not None and is_visible:
                _set_spine_position(ax_i.spines[side], ('outward', offset))

        # Set the ticks appropriately
        if bottom:
            ax_i.xaxis.tick_top()
        if top:
            ax_i.xaxis.tick_bottom()
        if left:
            ax_i.yaxis.tick_right()
        if right:
            ax_i.yaxis.tick_left()

        if trim:
            # clip off the parts of the spines that extend past major ticks
            xticks = ax_i.get_xticks()
            firsttick = np.compress(xticks >= ax_i.get_xlim()[0], xticks)[0]
            lasttick = np.compress(xticks <= ax_i.get_xlim()[-1], xticks)[-1]
            ax_i.spines['bottom'].set_bounds(firsttick, lasttick)
            ax_i.spines['top'].set_bounds(firsttick, lasttick)
            newticks = xticks.compress(xticks <= lasttick)
            newticks = newticks.compress(newticks >= firsttick)
            ax_i.set_xticks(newticks)

            yticks = ax_i.get_yticks()
            firsttick = np.compress(yticks >= ax_i.get_ylim()[0], yticks)[0]
            lasttick = np.compress(yticks <= ax_i.get_ylim()[-1], yticks)[-1]
            ax_i.spines['left'].set_bounds(firsttick, lasttick)
            ax_i.spines['right'].set_bounds(firsttick, lasttick)
            newticks = yticks.compress(yticks <= lasttick)
            newticks = newticks.compress(newticks >= firsttick)
            ax_i.set_yticks(newticks)


from astropy.io import fits
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm, PowerNorm
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_style('ticks')
# cmap = sns.cubehelix_palette(4, start=2.0, rot=0.2, dark=0.2, light=0.9, reverse=True)
cmap = sns.color_palette("cubehelix", 4)

# sns.set(context='paper')
def load_sdss_dr12(path):
    data_file = fits.open(path)
    # print(data_file[1].data.field)
    sdss_nam = data_file[1].data.field('SDSS_NAME')
    # f = 82900
    # print(sdss_nam[f:f+100])
    z_warning = data_file[1].data.field('zWarning')
    z = data_file[1].data.field('z_vi')
    # mask = np.logical_and(np.logical_and((z >= 1.0), (z <= 2.3)), (z_warning == 0))
    # mask =  np.ones(np.shape(z)).astype(bool)#(z_warning == 0)
    mask =  (z_warning == 0)
    z = z[mask]


    mi = data_file[1].data.field('MI')[mask]
    # print(len(mi))
    dgmi = data_file[1].data.field('DGMI')[mask]

    bands = ['u', 'g', 'r', 'i', 'z']
    data = {}
    length_data = []
    for i, k in enumerate(bands):
        data[k] = (data_file[1].data.field('PSFMAG')[:,i])[mask] - (data_file[1].data.field('EXTINCTION_RECAL')[:,i])[mask]
        length_data.append(len(data[k]))

    nam = np.array(['082045.38+130618.9', '115043.86-002354.1', '121940.36-010007.4', '123602.33-033129.9', '135425.24-001357.9', '143148.09+053558.0', '143748.28-014710.7'])
    # for n in nam:
    #     print(str(n))
    #     print([zip(i,k) for i,k in enumerate(sdss_nam) if str(n) == str(k)])
    # [print(k) for i,k in enumerate(sdss_nam)]
    u_obj = np.array([16.326057509461798, 16.857772100637099, 16.907661653092411, 16.805020474120248, 16.436534526124994, 16.058223981482136, 15.612695721407626])
    u_obj_sdss = np.array([16.28, 17.10, 17.22, 17.09, 16.98, 16.99, 15.86])
    g_obj = np.array([16.162734029157555, 16.942952211122567, 16.636137722929341, 16.743534687505424, 16.356290392576604, 16.880143790967175, 15.60611007862034])
    g_obj_sdss = np.array([16.13, 17.06, 16.92, 16.99, 16.78, 16.88, 15.76])
    r_obj = np.array([15.976060872654024, 16.911175035528821, 16.549549672624927, 16.595720735585481, 16.256037315738062, 16.765172373063898, 15.369589184790904])
    r_obj_sdss = np.array([15.91, 16.99, 16.82, 16.90, 16.67, 16.75, 15.48])
    i_obj = np.array([15.954628670601473, 16.647439017268361, 16.371386067198095, 16.364592262517355, 16.111299217989604, 16.585834842213593, 15.34670670905416])
    i_obj_sdss = np.array([15.87, 16.78, 16.63, 16.66, 16.51, 16.52, 15.41])
    z_obj = np.array([15.944986290774473, 16.501512015438244, 16.3521717641922, 16.358559488484381, 16.160066464330718, 16.350404442071188, 15.379652674965307])
    z_obj_sdss = np.array([15.83, 16.60, 16.61, 16.71, 16.51, 16.25, 15.41])
    mi_obj = np.array([-28.894596512440671, -29.49224989178208, -29.192578045151919, -29.530559286885733, -29.345333142602911, -29.759928110595922, -29.860980095511046])
    zz_obj = np.array([1.1242971107973012, 1.9798694976693576, 1.5830422701934881, 1.8463767030959246, 1.5123220764878522, 2.0997959346967061, 1.3089485173836708])


    print(np.mean(1 - u_obj / u_obj_sdss), np.std(1 - u_obj / u_obj_sdss))
    print(np.mean(1 - g_obj / g_obj_sdss), np.std(1 - g_obj / g_obj_sdss))
    print(np.mean(1 - r_obj / r_obj_sdss), np.std(1 - r_obj / r_obj_sdss))
    print(np.mean(1 - i_obj / i_obj_sdss), np.std(1 - i_obj / i_obj_sdss))
    print(np.mean(1 - z_obj / z_obj_sdss), np.std(1 - z_obj / z_obj_sdss))

    print(np.mean( u_obj -  u_obj_sdss), np.std( u_obj - u_obj_sdss))
    print(np.mean( g_obj -  g_obj_sdss), np.std( g_obj - g_obj_sdss))
    print(np.mean( r_obj -  r_obj_sdss), np.std( r_obj - r_obj_sdss))
    print(np.mean( i_obj -  i_obj_sdss), np.std( i_obj - i_obj_sdss))
    print(np.mean( z_obj -  z_obj_sdss), np.std( z_obj - z_obj_sdss))
    # pl.show()
    colors = ['z', 'g - i', 'i']
    gi = data['g'] - data['i']



    # data_color = np.array(zip(mi , gi))
    data_color = np.array(zip(z , gi, data['i']))




    # data_color = data_color[(np.logical_and(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5), data_color[:,0] < -28.0))]
    data_color = data_color[(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5))]
    # data = data[(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5))]
    data_color = data_color[(data_color[:,1] >= -5)]
    # data = data[(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5))]

    color_data = pd.DataFrame(data_color, columns=colors)


    # latexify()
    # Set up the subplot grid
    ratio = 5
    fig_width = 8

    golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_height = 1.0 * fig_width*golden_mean



    latexify(fig_width=fig_width, fig_height=fig_height)
    fig = pl.figure(figsize=(fig_width, fig_height))
    gs = pl.GridSpec(ratio + 1, ratio + 1)

    ax = fig.add_subplot(gs[1:, :-1])
    # ax_marg_x = fig.add_subplot(gs[0, :-1], sharex=ax)
    ax_marg_y = fig.add_subplot(gs[1:, -1], sharey=ax)

    x = np.array(color_data['z'])
    y = np.array(color_data['g - i'])
    imag = np.array(color_data['i'])


    import triangle

    # print(np.mean(x), np.std(y))
    print(np.mean(x), np.std(x))
    # p = sns.kdeplot(color_data, ax = ax, cmap=cmap, gridsize=10, linewidths = (0.5,))
    # p = sns.kdeplot(color_data, ax = ax, cmap=cmap, gridsize=500, linewidths = (0.5,), n_levels=25)

    # print(np.shape(data['i']))
    # print(np.shape(x))

    mask = (imag < 17.0) & ((1 < x) & (x < 2.3))
    # ax.scatter(mi_obj, g_obj - i_obj, s = 15, c=sns.xkcd_rgb["denim blue"], alpha = 0.7)
    triangle.hist2d(x, y, bins=200, ax=ax, smooth=0.3)
    ax.scatter(x[mask] , y[mask] ,  marker='o', s=10, facecolor=cmap[1], lw = 0, cmap=cmap, alpha= 1.0)
    ax.scatter(zz_obj, g_obj - i_obj, s = 25, c=cmap[2], alpha = 1.0)

    format_axes(ax)
    format_axes(ax_marg_y)
    pl.setp(ax_marg_y.get_yticklabels(), visible=False)
    pl.setp(ax_marg_y.yaxis.get_majorticklines(), visible=False)
    pl.setp(ax_marg_y.yaxis.get_minorticklines(), visible=False)
    pl.setp(ax_marg_y.xaxis.get_majorticklines(), visible=False)
    pl.setp(ax_marg_y.xaxis.get_minorticklines(), visible=False)
    pl.setp(ax_marg_y.get_xticklabels(), visible=False)
    ax_marg_y.xaxis.grid(False)
    despine(ax=ax_marg_y, bottom=True)
    sns.distplot(y, hist=False, kde=True, ax=ax_marg_y, kde_kws={"shade": True, "color": sns.xkcd_rgb["black"], "gridsize": 200, "alpha": 0.2},
                 vertical=True, axlabel=False)
    sns.distplot(g_obj - i_obj, hist=False, rug=True, kde=False, ax=ax_marg_y, rug_kws={"height": 1.5, "color": cmap[2], "alpha": 1.0},
                 vertical=True, axlabel=False)
    # sns.distplot(y[mask], hist=False, rug=True, kde=False, ax=ax_marg_y, rug_kws={"height": 0.5, "color": cmap[1], "alpha": 0.3},
    #              vertical=True, axlabel=False)
    sns.distplot(y[mask], hist=False, kde=True, ax=ax_marg_y, kde_kws={"shade": True, "color": cmap[1], "gridsize": 200, "alpha": 0.3},
                 vertical=True, axlabel=False)
    # sns.kdeplot(color_data, ax = ax, cmap='Greys', n_levels=50, norm=PowerNorm(gamma=0.3),
    #             shade=True, gridsize=100, linewidths = (0.5,), alpha=0.7)











    # ax.set_xlabel(r"M$_i$(z=2)")
    ax.set_xlabel(r"z")
    ax.set_ylabel(r"$g - i$")



    # ax.set_xlim((np.mean(x) - 5* np.std(x), np.mean(x) + 2* np.std(x)))
    ax.set_xlim((0.0, 3.5))
    ax.set_ylim((-0.5, 1.0))

    format_axes(ax)
    fig.tight_layout()
    pl.savefig('../documents/figs/color_comparison2.pdf')
    pl.show()





if __name__ == '__main__':
    # latexify()
    # path = '/Users/jselsing/nosync/sdss_quasar_catalog/DR12Q.fits'
    path = '/Users/jselsing/nosync/sdss_quasar_catalog/DR10Q_v2.fits'
    load_sdss_dr12(path=path)