#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')


import numpy as np
import glob
from scipy import interpolate
import matplotlib.pylab as pl
# use seaborn for nice default plot settings
import seaborn; seaborn.set_style('ticks')

from unred import ccm_unred,cardelli_reddening
from gen_methods import smooth,medfilt
from methods import common_wavelength

def main():
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'

    sdssobjects = glob.glob(root_dir+'*SDSS*/Telluric_corrected_science.dat')
    object_info_files = glob.glob(root_dir+'*SDSS*/Object_info.dat')


    obj_list =   [ 'SDSS0820+1306', 'SDSS1150-0023', 'SDSS1219-0100', 'SDSS1236-0331' , 'SDSS1354-0013',
                   'SDSS1431+0535', 'SDSS1437-0147']


    # obj_list = [ 'SDSS0820+1306']
    sdssobjects = [i for i in sdssobjects if i[-44:-31] in obj_list]
    object_info_files = [i for i in object_info_files if i[-29:-16] in obj_list]


    object_info = np.array([np.genfromtxt(i) for i in object_info_files])

    z_op = np.array([i[0] for i in object_info])
    z_sdss = np.array([i[1] for i in object_info])
    flag_z = np.array([i[2] for i in object_info])
    ebv = np.array([i[3] for i in object_info])

    redshifts = []
    for i,k in enumerate(flag_z):
        if k:
            redshifts.append(z_op[i])
        else:
            redshifts.append(z_sdss[i])


    # TODO Is there a better way to import the txt files?
    sdss_data_files = np.array([np.genfromtxt(i) for i in sdssobjects])

    wl = np.array([sdss_data_files[i][:,0] / (1 + redshifts[i]) for i in range(len(sdssobjects))])
    wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    flux = np.array([sdss_data_files[i][:,1] * (1 + redshifts[i]) for i in range(len(sdssobjects))])
    fluxerr = np.array([sdss_data_files[i][:,2] * (1 + redshifts[i]) for i in range(len(sdssobjects))])
    bp_map = np.array([sdss_data_files[i][:,3] for i in range(len(sdssobjects))])
    flux_cont = np.array([sdss_data_files[i][:,6] * (1 + redshifts[i]) for i in range(len(sdssobjects))])
    n_obj = len(obj_list)



    flux_johan = []
    fig, ax = pl.subplots(1, sharex=True)
    # for mask_ran in range(1500,9000, 1500):
    for mask_ran in [1420, 3000, 4000, 5600, 6800]:
    # for mask_ran in [6800]:
        # Construct common-size arrays:
        short = []
        tall = []
        for wl_i in wl:
            short.append(min(wl_i))
            tall.append(max(wl_i))
        short = min(short)
        tall = max(tall)
        step = 0.2 #CDELT
        wl_new = np.arange(short, tall, step)
        n_wl = len(wl_new)
        flux_new = np.zeros((n_obj,n_wl))
        flux_cont_new = np.zeros((n_obj,n_wl))
        fluxerr_new = np.zeros((n_obj,n_wl))
        bp_map_new = np.zeros((n_obj,n_wl))

        for n in range(n_obj):
            #de-reddening
            flux[n] = cardelli_reddening(wl_obs[n], flux[n], ebv[n])
            flux_cont[n] = cardelli_reddening(wl_obs[n], flux_cont[n], ebv[n])
            fluxerr[n] = cardelli_reddening(wl_obs[n], fluxerr[n], ebv[n])

        for n in range(n_obj):
            #Interpolate
            flux_new[n] = common_wavelength(wl[n], wl_new, flux[n])
            flux_cont_new[n] = common_wavelength(wl[n], wl_new, flux_cont[n])
            fluxerr_new[n] = common_wavelength(wl[n], wl_new, fluxerr[n], fill_value=1.0)
            bp_map_new[n] = common_wavelength(wl[n], wl_new, bp_map[n], fill_value=1.0)

        for n in range(n_obj):
            #Normalise
            mask = (wl_new > mask_ran) & (wl_new < mask_ran + 100)
            norm = np.median(flux_new[n][mask])
            flux_new[n] /= norm
            flux_cont_new[n] /= norm
            fluxerr_new[n] /= norm

        #Ensuring pixel usage blueward of Lya
        for i, k in enumerate(bp_map_new):
            mask_cont = (wl_new < 1216)
            (bp_map_new[i])[mask_cont] = 0

        #Saving ready spectra
        np.savetxt('test.dat', flux_new.transpose(), header="SDSS0820+1306 SDSS1150-0023 SDSS1219-0100 SDSS1236-0331 SDSS1354-0013 SDSS1431+0535 SDSS1437-0147")#, fmt = ['%5.1f', '%2.15E'] )



        #------------------------- Combination -------------------------
        # TODO Methodize this to avoid microchangning
        wmean = np.zeros(np.shape(wl_new))
        wmean_cont = np.zeros(np.shape(wl_new))
        mean = np.zeros(np.shape(wl_new))
        geo_mean = np.zeros(np.shape(wl_new))
        median = np.zeros(np.shape(wl_new))
        errofmean = np.zeros(np.shape(wl_new))
        errofwmean = np.zeros(np.shape(wl_new))
        std = np.zeros(np.shape(wl_new))
        CI_low = np.zeros(np.shape(mean))
        CI_high = np.zeros(np.shape(mean))
        for i, k in enumerate(flux_new.transpose()):
            mask = np.where(bp_map_new.transpose()[i] == 0)
            # mask = np.where(bp_map_new.transpose()[i] > -10000)
            # hej = True
            # if hej:
            if len(k[mask]) != 0:
                #Weighted mean
                weight = 1./(np.array(fluxerr_new.transpose()[i][mask])**2)
                wmean[i] = np.average(k[mask], axis = 0, weights = weight)
                wmean_cont[i] = np.average((flux_cont_new.transpose()[i])[mask], axis = 0, weights = weight)
                errofwmean[i] = np.sqrt(1./np.sum(np.array(fluxerr_new.transpose()[i][mask])**-2.,axis=0))

                #Mean
                mean[i] = np.mean(k[mask])
                errofmean[i] = np.sqrt(np.sum((fluxerr_new.transpose()[i][mask])**2))

                #Geometric mean
                from scipy.stats.mstats import gmean

                geo_mean[i] = gmean(k[mask])


                #Median
                median[i] = np.median(k[mask])
                # Bootstrapping to get confidence intervals
                # import scipy
                # import scikits.bootstrap as bootstrap
                # if len(k[mask]) > 1:
                #     CI_low[i], CI_high[i] = bootstrap.ci(data=k[mask], statfunction=np.median, n_samples=100, method='bca')
                # else:
                #     CI_low[i], CI_high[i] = np.median(k[mask]), np.median(k[mask])
                # if (i % 1000==0):
                #     print(i)


                std[i] = np.std(flux_new.transpose()[i][mask])
            elif len(k[mask]) == 0:

                mean[i] = 0
                errofmean[i] = 0
                wmean[i] = 0
                wmean_cont[i] = 0
                errofwmean[i] = 0
                geo_mean[i] = 0
                median[i] = 0
                CI_low[i], CI_high[i] = 0, 0
                std[i] = 0

        # ax.plot(wl_new, wmean_cont / np.median(wmean_cont), lw = 0.2)
        # ax[0].plot(wl_new, wmean / np.mean(wmean), lw = 0.5, label= str(mask_ran))
        # print(np.mean(wmean))
        # ii  = (wl_new > 6800) & (wl_new < 6800 + 100)
        # norm = np.median(wmean_cont[ii])
        # ax.plot(wl_new, wmean_cont / norm , lw = 0.5, label= str(mask_ran))

        ii  = (wl_new > 6800) & (wl_new < 6800 + 100)
        norm = np.median(geo_mean[ii])
        ax.plot(wl_new, geo_mean / norm , lw = 0.5, label= str(mask_ran))

        # ax.plot(wl_new, wmean_cont, lw = 0.5, label= str(mask_ran))

        # for n in range(n_obj):
            # plot individual spectra
            # ax.plot(wl_new, medfilt(flux_cont_new[n] , 51), lw=0.5, alpha = 0.5, color='grey')
        # pl.semilogy()
        # pl.show()


        # TODO Methodize this?
        # Calculating the number of spectra that goes into the composite and appending to n_spec.
        n_spec = np.zeros(np.shape(mean))
        spec = []
        for i, n in enumerate(flux_new.transpose()):
            mask = np.where(bp_map_new.transpose()[i] == 0)
            n_spec[i] = len((n[mask])[np.where((n[mask]) != 0)])
            if len(n[np.where(n != 0)]) == len(redshifts):
                spec.append(n / np.median(n))

        std = std / np.sqrt(n_spec)

        # ax[1].plot(wl_new, std, lw = 0.2)
        #Saving ready spectra

        #Saving to .dat file
        dt = [("wl", np.float64), ("wmean", np.float64) ]
        data = np.array(zip(wl_new, wmean), dtype=dt)
        file_name = "XSH-Composite_"+str(mask_ran)+""
        np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl wmean")#, fmt = ['%5.1f', '%2.15E'] )
        flux_johan.append(wmean.transpose())

    np.savetxt('test2.dat', zip(wl_new, *flux_johan), header="wl 1500 3000 4500 6000 7500")#, fmt = ['%5.1f', '%2.15E'] )

    ax.semilogy()
    ax.semilogx()
    # ax[1].semilogy()
    ax.legend()
    pl.show()

    # # Checking for normality
    # from matplotlib import pyplot as plt
    # import matplotlib.mlab as mlab
    # import scipy.stats as stats
    # p_val = []
    # # print(spec)
    # for i, k in enumerate(spec[10000:10110]):
    #     # pl.plot(n, '.' , hold=True)
    #     k = k / np.median(k)
    #     # print(k)
    #     k = np.hstack((k,np.mean(k)))
    #     #print(k)
    #     p_val.append((stats.normaltest(k)[1]))
    #     n, bins, patches = plt.hist(k, 10, hold=True)
    #     mu = np.mean(k)
    #     sigma = np.std(k)
    #     plt.plot(bins, mlab.normpdf(bins, mu, sigma), hold=True)
    # print(np.mean(p_val))
    # pl.xlabel(r'Normalised value [input/median(input)]')
    # pl.ylabel(r'Arbitrary scale')
    # pl.savefig("../documents/normality.pdf", dpi= 150)
    # pl.show()




    #Saving to .dat file
    dt = [("wl", np.float64), ("mean", np.float64), ("mean_error", np.float64), ("wmean", np.float64),
          ("wmean_error", np.float64), ("geo_mean", np.float64), ("median", np.float64), ("n_spec", np.float64),
          ("std", np.float64), ("wmean_cont", np.float64) ]
    data = np.array(zip(wl_new, mean, errofmean, wmean, errofwmean, geo_mean, median, n_spec, std, wmean_cont ), dtype=dt)
    file_name = "Composite"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl mean mean_error wmean wmean_error geo_mean median n_spec std wmean_cont")#, fmt = ['%5.1f', '%2.15E'] )

    # #Saving to .dat file
    # dt = [("wl", np.float64), ("wmean", np.float64) ]
    # data = np.array(zip(wl_new, wmean), dtype=dt)
    # file_name = "XSH-Composite_7500"
    # np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl wmean")#, fmt = ['%5.1f', '%2.15E'] )

    #Saving to .dat file
    dt = [("wl", np.float64), ("wmean_cont", np.float64) ]
    data = np.array(zip(wl_new, wmean_cont), dtype=dt)
    file_name = "XSH-Composite_cont"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl wmean_cont")#, fmt = ['%5.1f', '%2.15E'] )

if __name__ == '__main__':
    main()