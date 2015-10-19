#!/usr/bin/env python
# -*- coding: utf-8 -*-



from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


from methods import latexify, format_axes, gauss


import numpy as np

if __name__ == '__main__':
    dat = np.genfromtxt('data/regularised.dat')
    print(np.shape(dat))
    n_test = (dat[5800 : 5900,:])
    # Checking for normality
    from matplotlib import pyplot as plt
    import seaborn as sns; sns.set_style('ticks')
    cmap = sns.color_palette("cubehelix", 6)
    import scipy.stats as stats
    import statsmodels.api as sm
    p_val = []

    #Plotting
    ratio = (1.0 + np.sqrt(5.0))/2.0
    latexify(columns=2)
    fig, ax = plt.subplots()
    for i, k in enumerate(n_test):

        p_val.append((stats.shapiro(k)[1]))

    mtest = np.mean(n_test, axis = 0)
    print(mtest)

    sm.qqplot(mtest, fit=True, line='45', ax=ax)

    print(np.mean(p_val))

    format_axes(ax)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    fig.tight_layout()
    plt.savefig("../documents/figs/normality.pdf", dpi= 150)
    plt.show(block=True)