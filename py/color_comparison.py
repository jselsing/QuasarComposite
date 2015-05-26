#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"

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


import matplotlib

from math import sqrt
SPINE_COLOR = 'gray'

def latexify(fig_width=None, fig_height=None, columns=1):
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
              'axes.labelsize': 11, # fontsize for x and y labels (was 10)
              'axes.titlesize': 11,
              'font.size': 11, # was 10
              'legend.fontsize': 11, # was 10
              'xtick.labelsize': 11,
              'ytick.labelsize': 11,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }




    matplotlib.rcParams.update(params)


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

from astropy.io import fits
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm, PowerNorm
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_style('ticks')
cmap = sns.cubehelix_palette(10, start=2, rot=0, dark=0.3, light=0.9, reverse=True, as_cmap=True)
# sns.set(context='paper')
def load_sdss_dr12(path):
    data_file = fits.open(path)
    # print(data_file[1].data.field)
    z_warning = data_file[1].data.field('zWarning')
    z = data_file[1].data.field('z_vi')
    # mask = np.logical_and(np.logical_and((z >= 1.0), (z <= 2.3)), (z_warning == 0))
    mask =  (z_warning == 0)
    # print(z[mask])
    z = z[mask]
    print(len(z))

    mi = data_file[1].data.field('MI')[mask]
    print(len(mi))
    dgmi = data_file[1].data.field('DGMI')[mask]

    bands = ['u', 'g', 'r', 'i', 'z']
    data = {}
    length_data = []
    for i, k in enumerate(bands):
        data[k] = (data_file[1].data.field('PSFMAG')[:,i])[mask] - (data_file[1].data.field('EXTINCTION_RECAL')[:,i])[mask]
        length_data.append(len(data[k]))


    u_obj = np.array([16.776679059420339, 17.695431195769537, 17.375712062291065, 17.491023145887873, 17.159539776215617, 17.584843370477948, 16.257372339772438])
    g_obj = np.array([16.7175770449848, 17.66287889255684, 17.264278555286204, 17.423701102087442, 17.168422183077446, 17.590939752133799, 16.155141214036895])
    r_obj = np.array([16.581391280385212, 17.420886487129813, 17.164700226920623, 17.309943991447987, 17.038907610049456, 17.338833016940775, 15.915406330517975])
    i_obj = np.array([16.846063280799477, 17.766023210621618, 17.420903544404176, 17.600066414834963, 17.366033651901112, 18.016017836357896, 16.200013990848682])
    z_obj = np.array([16.795833639704391, 21.065331593637801, 17.507279309060571, 19.021921387107618, 17.512653633591903, 25.971125223599778, 16.175981506541454])
    mi_obj = np.array([-28.894596512440671, -29.49224989178208, -29.192578045151919, -29.530559286885733, -29.345333142602911, -29.759928110595922, -29.860980095511046])
    zz_obj = np.array([1.1242971107973012, 1.9798694976693576, 1.5830422701934881, 1.8463767030959246, 1.5123220764878522, 2.0997959346967061, 1.3089485173836708])
    # pl.show()
    colors = ['M_i', 'g - i']
    gi = data['g'] - data['i']



    # data_color = np.array(zip(mi , gi))
    data_color = np.array(zip(z , gi))



    # data_color = data_color[(np.logical_and(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5), data_color[:,0] < -28.0))]
    data_color = data_color[(np.logical_and(data_color[:,0] > -40.0, data_color[:,1] >= -5))]
    data_color = data_color[(data_color[:,1] >= -5)]


    color_data = pd.DataFrame(data_color, columns=colors)


    latexify()
    # Set up the subplot grid
    ratio = 5
    fig_width = 6

    golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_height = fig_width*golden_mean


    fig = pl.figure(figsize=(fig_width, fig_height))
    gs = pl.GridSpec(ratio + 1, ratio + 1)

    ax = fig.add_subplot(gs[1:, :-1])
    # ax_marg_x = fig.add_subplot(gs[0, :-1], sharex=ax)
    ax_marg_y = fig.add_subplot(gs[1:, -1], sharey=ax)

    x = np.array(color_data['M_i'])
    y = np.array(color_data['g - i'])


    import triangle

    # print(np.mean(x), np.std(y))
    print(np.mean(x), np.std(x))
    # p = sns.kdeplot(color_data, ax = ax, cmap=cmap, gridsize=10, linewidths = (0.5,))
    # p = sns.kdeplot(color_data, ax = ax, cmap=cmap, gridsize=500, linewidths = (0.5,), n_levels=25)


    mask = (data['i'] < 17.0)
    # ax.scatter(mi_obj, g_obj - i_obj, s = 15, c=sns.xkcd_rgb["denim blue"], alpha = 0.7)
    triangle.hist2d(x, y, bins=200, ax=ax, smooth=0.3)
    ax.scatter(x[mask] , y[mask] ,  marker='o', s=10, facecolor=sns.xkcd_rgb["steel grey"], lw = 0, cmap=cmap, alpha= 1.0)
    ax.scatter(zz_obj, g_obj - i_obj, s = 25, c=sns.xkcd_rgb["denim blue"], alpha = 1.0)

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
    sns.distplot(y, hist=False, kde=True, ax=ax_marg_y, kde_kws={"shade": True, "color": sns.xkcd_rgb["black"], "gridsize": 200},
                 vertical=True, axlabel=False)
    sns.distplot(g_obj - i_obj, hist=False, rug=True, kde=False, ax=ax_marg_y, rug_kws={"height": 1.5, "color": sns.xkcd_rgb["denim blue"], "alpha": 1.0},
                 vertical=True, axlabel=False)
    sns.distplot(y[mask], hist=False, rug=True, kde=False, ax=ax_marg_y, rug_kws={"height": 0.5, "color": sns.xkcd_rgb["steel grey"], "alpha": 0.3},
                 vertical=True, axlabel=False)

    # sns.kdeplot(color_data, ax = ax, cmap='Greys', n_levels=50, norm=PowerNorm(gamma=0.3),
    #             shade=True, gridsize=100, linewidths = (0.5,), alpha=0.7)











    # ax.set_xlabel(r"M$_i$(z=2)")
    ax.set_xlabel(r"z")
    ax.set_ylabel("g - i")



    # ax.set_xlim((np.mean(x) - 5* np.std(x), np.mean(x) + 2* np.std(x)))
    ax.set_xlim((0.0, 3.5))
    ax.set_ylim((-0.5, 1.0))


    fig.tight_layout()
    pl.savefig('../documents/figs/color_comparison2.pdf', dpi=2400)
    pl.show()





if __name__ == '__main__':
    latexify()
    # path = '/Users/jselsing/nosync/sdss_quasar_catalog/DR12Q.fits'
    path = '/Users/jselsing/nosync/sdss_quasar_catalog/DR10Q_v2.fits'
    load_sdss_dr12(path=path)