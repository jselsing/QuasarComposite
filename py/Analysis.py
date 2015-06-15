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

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


#
#

def main():
    from matplotlib import rc_file
    rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
    import numpy as np

    import matplotlib.pylab as pl
    # use seaborn for nice default plot settings
    import seaborn; seaborn.set_style('ticks')


    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    #
    # data_file = np.genfromtxt(root_dir+'untellcorrComposite.dat')
    #
    #
    # wluntell = data_file[:,0]
    # meanuntell = data_file[:,1]
    # err_meanuntell = data_file[:,2]
    # wmeanuntell = data_file[:,3]
    # err_wmeanuntell = data_file[:,4]
    # geo_meanuntell = data_file[:,5]
    # medianuntell = data_file[:,6]
    # n_specuntell = data_file[:,7]
    # stduntell = data_file[:,8]
    # # wl_sdss = data_file[:,9]
    # # mean_sdss = data_file[:,10]
    # wmean_contuntell = data_file[:,9]






    data_file = np.genfromtxt(root_dir+'Composite.dat')
    # data_file = np.genfromtxt(root_dir+'untellcorrComposite.dat')


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
    # par_guess3 = [1, -1.0, -0.00001]

    # pl.plot(wl, power_law2(wl, par_guess))
    # pl.show()
    wmean[np.where(np.isnan(wmean) == True)] = 0

    mask = (wl > 1300) & (wl < 1350) | (wl > 1425) & (wl < 1475) | (wl > 5500) & (wl < 5800) | (wl > 7300) & (wl < 7500)
    # mask = (wl > 1080) & (wl < 1150) | (wl > 1350) & (wl < 1365) | (wl > 4200) & (wl < 4230) | (wl > 5500) & (wl < 6035) | (wl > 7800) & (wl < 7950)


    # mask = np.where(wl > 0)
    # print(mask)
    popt, pcov = optimize.curve_fit(power_law, wl[mask], wmean_cont[mask], p0=par_guess, sigma=err_wmean[mask], absolute_sigma=True)
    popt2, pcov2 = optimize.curve_fit(power_law2, wl[mask], wmean_cont[mask], p0=par_guess2, sigma=err_wmean[mask], absolute_sigma=True)
    # popt3, pcov3 = optimize.curve_fit(power_law3, wl[mask], wmean[mask], p0=par_guess3)


    print(*popt)
    print(*popt2)
    # print(*popt3)


    # Smoothing for presentation purposes
    # from gen_methods import smooth,medfilt
    # from xshoo.binning import binning1d
    # bins = 5
    # wl = binning1d(wl, bins)
    # mean, std = binning1d(mean, bins, err=errofmean)
    # mean = medfilt(mean, 5)
    # errofmean = medfilt(errofmean, 5)
    # wmean = medfilt(wmean, 35)
    # errofwmean = medfilt(errofwmean, 5)
    # geo_mean = medfilt(geo_mean, 5)
    # median = medfilt(median, 5)
    # median = smooth(np.array(median), window_len=5, window='hanning')
    # wmean = smooth(np.array(wmean), window_len=5, window='hanning')

    #
    #
    #
    # fig, (ax1, ax2, ax3, ax4) = pl.subplots(ncols=1, nrows=4, sharex=True)
    # linewidth = 0.1
    #
    # ax4.plot(wl,medfilt( n_spec, 151), lw = linewidth, color = 'black', label='Number of spectra')
    # ax4.semilogx()
    # ax4.set_xlim((1000, 11000 ))
    # # pl.xlabel(r'Rest wavelength [\AA]')
    # # pl.ylabel(r'Number of spectra')
    # # pl.xlim((1000, 11000 ))
    # # ymin = 0# 0.3
    # # ymax = 10#155
    # # ax4.set_ylim = ((ymin, ymax))
    # # pl.savefig("../documents/figs/number_spec.pdf", dpi= 150)
    # # pl.show()
    # # pl.clf()
    # #
    # #
    # #
    # ax2.plot(wl, medfilt(wmean/err_wmean, 5), lw=linewidth, color = 'black' , linestyle = 'steps-mid', label = 'Signal-to-noise')
    # ax2.semilogx()
    # # pl.xlabel(r'Rest wavelength [\AA]')
    # # pl.ylabel(r'S/N')
    # ax2.set_xlim((1000, 11000 ))
    # # pl.savefig("../documents/figs/signal_to_noise.pdf", dpi= 150)
    # # pl.show()
    # # pl.clf()
    #
    # ax3.plot(wl, medfilt(std, 15), lw=linewidth, color = 'black' , linestyle = 'steps-mid', label='variance')
    # # pl.xlabel(r'Rest wavelength [\AA]')
    # # pl.ylabel(r'Standard deviation of constituent spectra')
    # ax3.semilogy()
    # ax3.semilogx()
    # ax3.set_xlim((1000, 11000 ))
    # ax3.set_ylim((1e-3, 1e2))
    # # pl.savefig("../documents/figs/std.pdf", dpi= 150)
    # # pl.show()
    # # pl.clf()
    #
    #
    #
    # #Plotting
    #
    # ax1.plot(wl, power_law2(wl, *popt2) , 'b--')
    # ax1.plot(wl, power_law(wl, *popt) , 'b--')
    #
    # # pl.plot(wl, wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # ax1.plot(wl, wmean_cont, color = 'black', lw = linewidth, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # ax1.semilogy()
    # ax1.semilogx()
    # # pl.plot(wl, mean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter mean composite')
    # # pl.plot(wl_sdss, mean_sdss, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # # pl.plot(wl, geo_mean, color = 'red', lw = 0.5, linestyle = 'steps-mid', label='X-shooter geometric composite')
    # # pl.plot(wl, median, color = 'green', lw = 0.5, linestyle = 'steps-mid', label='X-shooter median composite')
    # # pl.plot(wl, err_mean, color = 'black', lw = 0.5, linestyle = 'steps-mid')
    # # pl.plot(wl, wmean / err_wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid')
    #
    #
    #
    #
    #
    #
    # #Overplot lines
    # fit_line_positions = np.genfromtxt('plotlinelist.txt', dtype=None)
    # linelist = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    # linelist = np.array(linelist)
    #
    # ymin = 0.3# 0.3
    # ymax = 150#155
    # ax1.set_ylim((ymin, ymax))
    # for p in range(len(fit_line_positions)):
    #     xcoord = linelist[p]
    #     mask = (wl > xcoord - 1) & (wl < xcoord + 1)
    #     y_val = np.mean(wmean_cont[mask])
    #     ax1.axvline(x=xcoord, ymin=(np.log(y_val) + 1.30) / (( np.log(ymax) - np.log(ymin))), color='red',linestyle='dashed', lw=0.75)
    #     trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
    #     ann = ax1.annotate(fit_line_positions[p,][0], xy=(xcoord, 1.05 ), fontsize='x-small', rotation = 90, xycoords=trans)
    #
    #
    # # ------------------------- Auxilliary products -------------------------
    # # Johans Composite:
    # vandenberk = np.loadtxt('Glikman.data')
    # mask = (wl > 3600) & (wl < 3700)
    # maskVB = (vandenberk[:,0] > 3600) & (vandenberk[:,0] < 3700)
    # norm = np.median(vandenberk[:,1][maskVB]) / np.median(mean[mask])
    # # Vanden Berk:
    # vandenberk2 = np.loadtxt('composite.txt')
    # mask = (wl > 3600) & (wl < 3700)
    # maskVB = (vandenberk2[:,0] > 3600) & (vandenberk2[:,0] < 3700)
    # norm2 = np.median(vandenberk2[:,1][maskVB]) / np.median(mean[mask])
    # # pl.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Glikman Composite')
    # # pl.plot(vandenberk2[:,0],vandenberk2[:,1]/norm2,drawstyle='steps-mid',label='Vanden Berk Composite')
    #
    #
    # # pl.xlim((1150,10500))
    # # pl.ylim((ymin,ymax))
    #
    # # pl.xlim((0.95*min(wl), 1.01*max(wl)))
    # pl.xlabel(r'Rest Wavelength  [\r{A]')
    # ax1.set_ylabel(r'Normalised Flux')
    # ax2.set_ylabel(r'Ratio')
    # ax3.set_ylabel(r'Variance')
    # ax4.set_ylabel(r'Number of spectra')
    #
    #
    #
    #
    # pl.xlim((1000, 11000 ))
    # pl.semilogx()
    # pl.loglog()
    # ax1.legend()
    # ax2.legend()
    # ax3.legend()
    # ax4.legend()
    # pl.tight_layout()
    # fig.subplots_adjust(hspace=0)
    # fig = pl.gcf()
    # fig.set_size_inches(10,10)
    # fig.savefig("../documents/figs/Combined.pdf", dpi= 150)
    # pl.show()
    #
    #
    # #Plotting
    # # pl.plot(wl, power_law3(wl, *popt3) , 'b--')
    # pl.plot(wl, power_law2(wl, *popt2) , 'b--')
    # pl.plot(wl, power_law(wl, *popt) , 'b--')
    #
    # # pl.plot(wl, wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # pl.plot(wl, wmean_cont, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # # pl.plot(wl, mean, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter mean composite')
    # # pl.plot(wl_sdss, mean_sdss, color = 'black', lw = 0.5, linestyle = 'steps-mid', label='X-shooter wmean composite')
    # # pl.plot(wl, geo_mean, color = 'red', lw = 0.5, linestyle = 'steps-mid', label='X-shooter geometric composite')
    # # pl.plot(wl, median, color = 'green', lw = 0.5, linestyle = 'steps-mid', label='X-shooter median composite')
    # # pl.plot(wl, err_mean, color = 'black', lw = 0.5, linestyle = 'steps-mid')
    # # pl.plot(wl, wmean / err_wmean, color = 'black', lw = 0.5, linestyle = 'steps-mid')
    #
    #
    #


    # #Overplot lines
    # fit_line_positions = np.genfromtxt('plotlinelist.txt', dtype=None)
    # linelist = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    # linelist = np.array(linelist)
    #
    # ymin = 0.3# 0.3
    # ymax = 150#155
    # # ax1.set_ylim((ymin, ymax))
    # for p in range(len(fit_line_positions)):
    #     xcoord = linelist[p]
    #     mask = (wl > xcoord - 1) & (wl < xcoord + 1)
    #     y_val = np.mean(wmean_cont[mask])
    #     pl.axvline(x=xcoord, ymin=(np.log(y_val) + 1.30) / (( np.log(ymax) - np.log(ymin))), color='red',linestyle='dashed', lw=0.75)
    #     trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
    #     ann = pl.annotate(fit_line_positions[p,][0], xy=(xcoord, 1.05 ), fontsize='x-small', rotation = 90, xycoords=trans)
    #     # ax1.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 2.0 ),fontsize='x-small', rotation = 90)
    #
    #
    # # ------------------------- Auxilliary products -------------------------
    # # Johans Composite:
    # vandenberk = np.loadtxt('Glikman.data')
    # mask = (wl > 3600) & (wl < 3700)
    # maskVB = (vandenberk[:,0] > 3600) & (vandenberk[:,0] < 3700)
    # norm = np.median(vandenberk[:,1][maskVB]) / np.median(mean[mask])
    # # Vanden Berk:
    # vandenberk2 = np.loadtxt('composite.txt')
    # mask = (wl > 3600) & (wl < 3700)
    # maskVB = (vandenberk2[:,0] > 3600) & (vandenberk2[:,0] < 3700)
    # norm2 = np.median(vandenberk2[:,1][maskVB]) / np.median(mean[mask])
    # # pl.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Glikman Composite')
    # # pl.plot(vandenberk2[:,0],vandenberk2[:,1]/norm2,drawstyle='steps-mid',label='Vanden Berk Composite')
    #
    #
    #
    # pl.xlabel(r'Rest Wavelength  [\r{A}]')
    # pl.ylabel(r'Normalised Flux')
    #
    #
    #
    #
    # pl.xlim((1000, 11000 ))
    # # pl.semilogx()
    # pl.loglog()
    # # ax1.legend()
    # # ax2.legend()
    # # ax3.legend()
    # # ax4.legend()
    # # pl.tight_layout()
    # # fig.subplots_adjuRst(hspace=0)
    # fig = pl.gcf()
    # # fig.set_size_inches(3,6)
    # fig.savefig("../documents/figs/Composite.pdf", dpi= 150)
    # pl.show()

    from gen_methods import medfilt, smooth
    def hist(rawData,xRange,nBins=10,mode='lin'):

        """histogram using linear binning of supplied data

        Input:
        rawData	-- list containing data to be binned
        xRange  -- lower(incl)/upper(excl) boundary for numerical values
        nBin    -- desired number of bins (default =10)
        mode	-- binning type (possible choices: lin, log)

        Returns: (nothing)
        """
        from math   import sqrt,floor,log,exp
        h = [0]*nBins
        xMin=float(xRange[0])
        xMax=float(xRange[1])

        if mode == 'lin':
            dx = (xMax-xMin)/nBins
            def binId(val):   return int(floor((val-xMin)/dx))
            def bdry(bin):	  return xMin+bin*dx, xMin+(bin+1)*dx
            def GErr(q,n,dx): return sqrt(q*(1-q)/(N-1))/dx

        elif mode == 'log':
            dx = log(xMax/xMin)/nBins
            def binId(val):   return int(floor(log(val/xMin)/dx))
            def bdry(bin):	  return xMin*exp(bin*dx), xMin*exp((bin+1)*dx)
            def GErr(q,n,dx): return "##"

        for value in rawData:
            if 0<=binId(value)<nBins:
              h[binId(value)] += 1

        N=sum(h)
        binned = []
        for bin in range(nBins):
            hRel   = float(h[bin])/N
            low,up = bdry(bin)
            binned.append(low)
            width  = up-low
            # print(low, up, hRel/width, GErr(hRel,N,width))
        return binned

    nbins = len(wl)
    log_binned_wl = hist(wl,[min(wl),max(wl)], nbins,'log')
    from scipy.interpolate import InterpolatedUnivariateSpline
    sps = InterpolatedUnivariateSpline(wl, wmean_cont)




    # import matplotlib.pylab as pl
    # import lineid_plot
    # # use seaborn for nice default plot settings
    # import seaborn; seaborn.set_style('ticks')
    # # cmap = seaborn.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True)
    # seaborn.set_palette('muted')
    # # deep, muted, bright, pastel, dark, colorblind
    # fig, ax = pl.subplots(1 , figsize=(12, 4))
    # ax.plot(log_binned_wl, medfilt(sps(log_binned_wl) , 21),
    #         lw = 1.0, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')
    # # ax.plot(wl, wmean_cont,
    # #         lw = 0.5, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')
    #
    # ax.plot(wl, power_law(wl, *popt),
    #         linestyle='dashed', label ='Power law fit')
    # ax.plot(wl, power_law2(wl, *popt2),
    #         linestyle='dashed', label ='Power law fit')
    #
    #
    #
    #
    # sdss_compo = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/sdss_compo.dat')
    # sdss_wl = sdss_compo[:,0]
    # sdss_flux = sdss_compo[:, 1]
    #
    # norm_reg = 1430
    #
    # mask = (wl > norm_reg) & (wl < norm_reg + 20)
    # norm1 = np.median(wmean_cont[mask])
    # # mask = (log_binned_wl > norm_reg) & (log_binned_wl < norm_reg + 10)
    # # norm1 = np.median(medfilt(sps(log_binned_wl) , 21)[mask])
    # print(norm1)
    # mask = (sdss_wl > norm_reg) & (sdss_wl < norm_reg + 20)
    # norm2 = np.median(sdss_flux[mask])
    # print(norm2)
    #
    #
    # ax.plot(sdss_wl, sdss_flux * (norm1/norm2),
    #         linestyle='solid', label ='SDSS composite')
    # #Overplot lines
    # fit_line_positions = np.genfromtxt('plotlinelist.txt', dtype=None)
    #
    # linelist = []
    # linenames = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    #     linenames.append(n[0])
    # # linelist = np.array(linelist)
    #
    # ax.set_xlim((950, 12500))
    # ax.set_ylim((0.09, 25))
    # pl.semilogy()
    # # pl.semilogx()
    #
    #
    # # Formatting axes
    # import matplotlib as mpl
    # ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_xticks([1000, 2000, 5000, 10000])
    # ax.get_xaxis().tick_bottom()
    # # ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    #
    # ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_yticks([0.1, 1, 10, 20])
    #
    # lineid_plot.plot_line_ids(wl, wmean_cont, linelist, linenames, ax=ax)
    # for i in ax.lines:
    #     if '$' in i.get_label():
    #         i.set_alpha(0.3)
    # # pl.legend()
    # pl.savefig('../documents/figs/compo.pdf')
    # pl.show()






    # n_spec[n_spec == 0] = 1
    # wmean_cont[wmean_cont == 0] = 1
    # # std[std < 0.005] = 0.005
    # std[n_spec == 1] = 0
    # std[std < 0.005] = 0.005
    #
    # # print(std/wmean_cont)
    # sps = InterpolatedUnivariateSpline(wl, std)#/wmean_cont)
    #
    #
    # import matplotlib.pylab as pl
    # import lineid_plot
    # # use seaborn for nice default plot settings
    # import seaborn; seaborn.set_style('ticks')
    # # cmap = seaborn.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True)
    # seaborn.set_palette('muted')
    # # deep, muted, bright, pastel, dark, colorblind
    # fig, ax = pl.subplots(1 , figsize=(12, 4))
    # ax.plot(log_binned_wl, medfilt(sps(log_binned_wl) , 11),
    #         lw = 0.5, alpha=1.0, linestyle = 'steps-mid', label='Composite variability')
    #
    #
    # #Overplot lines
    # fit_line_positions = np.genfromtxt('plotlinelist.txt', dtype=None)
    #
    # linelist = []
    # linenames = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    #     linenames.append(n[0])
    # # linelist = np.array(linelist)
    #
    # ax.set_xlim((950, 11500))
    # # ax.set_ylim((0.001, 25))
    # ax.set_ylim((0.001, 1))
    # # pl.semilogy()
    # pl.semilogx()
    #
    #
    # # Formatting axes
    # import matplotlib as mpl
    # ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_xticks([1000, 2000, 5000, 10000])
    # ax.get_xaxis().tick_bottom()
    # # ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    #
    # ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    # # ax.set_yticks([0.1, 1, 10, 20])
    #
    # # lineid_plot.plot_line_ids(wl, std/wmean_cont, linelist, linenames, ax=ax)
    # lineid_plot.plot_line_ids(wl, std, linelist, linenames, ax=ax)
    # for i in ax.lines:
    #     if '$' in i.get_label():
    #         i.set_alpha(0.3)
    # # pl.legend()
    # pl.savefig('../documents/figs/var.pdf')
    # pl.show()





    # out = std#/wmean_cont
    #
    # #Saving to .dat file
    # dt = [("wl", np.float64), ("out", np.float64) ]
    # data = np.array(zip(wl, out), dtype=dt)
    # file_name = "Spectral variability"
    # np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl variability")#, fmt = ['%5.1f', '%2.15E'] )



    import matplotlib.pylab as pl
    import lineid_plot
    # use seaborn for nice default plot settings
    import seaborn; seaborn.set_style('ticks')
    # cmap = seaborn.cubehelix_palette(3, start=2, rot=0.0, dark=0.2, light=.8, reverse=True)
    # seaborn.set_palette(cmap)

    # deep, muted, bright, pastel, dark, colorblind

    fig, ax = pl.subplots(1 , figsize=(12, 4))
    # print(ax)
    hej = medfilt(sps(log_binned_wl) , 1)
    print(hej)
    ax.plot(log_binned_wl, hej,
            lw = 1.0, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')

    # hejn = hej[3000]
    # sps = InterpolatedUnivariateSpline(wl, wmeanuntell)
    # hej2 = medfilt(sps(log_binned_wl) , 21)
    # hej2n = hej2[3000]
    # norm = ( hejn/ hej2n)
    # ax.plot(log_binned_wl, hej2 * norm ,
    #         lw = 1.0, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')
    #
    # mask = (log_binned_wl > 1300) & (log_binned_wl < 5000)
    # print(np.mean((hej / hej2 * norm)[mask]), np.std((hej / hej2 * norm)[mask]))
    # ax.plot(log_binned_wl[mask], (hej / hej2 * norm)[mask] ,
    #         lw = 1.0, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')


    # ax.plot(wl, wmean_cont,
    #         lw = 0.5, alpha=1.0, linestyle = 'steps-mid', label='X-shooter mean composite')

    ax.plot(wl, power_law(wl, *popt),
            linestyle='dashed', label ='Power law fit')
    ax.plot(wl, power_law2(wl, *popt2),
            linestyle='dashed', label ='Power law fit')




    sdss_compo = np.genfromtxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/sdss_compo.dat')
    sdss_wl = sdss_compo[:,0]
    sdss_flux = sdss_compo[:, 1]

    norm_reg = 1430

    mask = (wl > norm_reg) & (wl < norm_reg + 20)
    norm1 = np.median(wmean_cont[mask])
    # mask = (log_binned_wl > norm_reg) & (log_binned_wl < norm_reg + 10)
    # norm1 = np.median(medfilt(sps(log_binned_wl) , 21)[mask])
    print(norm1)
    mask = (sdss_wl > norm_reg) & (sdss_wl < norm_reg + 20)
    norm2 = np.median(sdss_flux[mask])
    print(norm2)


    ax.plot(sdss_wl, sdss_flux * (norm1/norm2),
            linestyle='solid', label ='Full sample SDSS composite')
    #Overplot lines
    fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)

    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(n[1])
        linenames.append(n[0])
    # linelist = np.array(linelist)


    pl.semilogy()
    # pl.semilogx()

    # Formatting axes
    import matplotlib as mpl
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xticks([1100, 3000, 5000])
    ax.get_xaxis().tick_bottom()
    # ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50])

    pl.legend(loc=3)
    ax.set_xlim((1100, 11000))
    ax.set_ylim((0.1, 50))




    lineid_plot.plot_line_ids(wl, wmean_cont, linelist, linenames, ax=ax)
    for i in ax.lines:
        if '$' in i.get_label():
            i.set_alpha(0.5)

    # pl.savefig('../documents/figs/compo_full_sample.pdf')
    pl.show()






if __name__ == '__main__':
    main()