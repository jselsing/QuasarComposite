#!/usr/bin/env python
# -*- coding: utf-8 -*-



from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')

from methods import latexify, format_axes, gauss


import numpy as np

if __name__ == '__main__':
    dat = np.genfromtxt('data/regularised.dat')
    print(np.shape(dat))
    n_test = (dat[5800 : 5900,:-1])
    # Checking for normality
    from matplotlib import pyplot as plt
    import seaborn as sns; sns.set_style('ticks')
    cmap = sns.color_palette("cubehelix", 6)
    import scipy.stats as stats
    import statsmodels.api as sm
    p_val = []

    #Plotting
    ratio = (1.0 + np.sqrt(5.0))/2.0
    latexify(fig_width=5*ratio, fig_height=5)
    fig, ax = plt.subplots(figsize=(5*ratio, 5))
    for i, k in enumerate(n_test):
        # k = np.hstack((k, [np.median(k)]))
        # pl.plot(n, '.' , hold=True)
        # k = k / np.median(k)
        # print(k)
        # k = np.hstack((k,np.mean(k)))
        # print(k)
        # print(k)
        p_val.append((stats.shapiro(k)[1]))
    # n, bins, patches = plt.hist(k, 10, hold=True)
        # mu = np.mean(k)
        # sigma = np.std(k)
        # plt.plot(bins, mlab.normpdf(bins, mu, sigma), hold=True)
    mtest = np.mean(n_test, axis = 0)
    # mtest = n_test
    # print(k)
    # stats.probplot(k, dist="norm", plot=plt)
    sm.qqplot(mtest, fit=True, line='45', ax=ax)
    print(np.mean(p_val))
    # latexify(fig_width=5*ratio, fig_height=5)
    format_axes(ax)

    # pl.xlabel(r'Normalised value [input/median(input)]')
    # pl.ylabel(r'Arbitrary scale')
    fig.tight_layout()
    plt.savefig("../documents/figs/normality.pdf", dpi= 150)
    plt.show(block=True)