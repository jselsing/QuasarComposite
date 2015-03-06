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
from unred import ccm_unred,cardelli_reddening


def main():
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'

    sdssobjects = glob.glob(root_dir+'*SDSS*/Telluric_corrected_science.dat')
    object_info_files = glob.glob(root_dir+'*SDSS*/Object_info.dat')


    obj_list =   [ 'SDSS0820+1306', 'SDSS1150-0023', 'SDSS1219-0100', 'SDSS1236-0331' , 'SDSS1354-0013',
                   'SDSS1431+0535', 'SDSS1437-0147']
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
    # redshifts = [1.1250, 1.9798, 1.5826, 1.8458, 1.5123, 2.0998, 1.3092]
    # ebv = [0.0253, 0.0201, 0.0279, 0.0268, 0.0330, 0.0274, 0.0378]
    # TODO Is there a better way to import the txt files?
    sdss_data_files = np.array([np.genfromtxt(i) for i in sdssobjects])

    wl = np.array([sdss_data_files[i][:,0] / (1 + redshifts[i]) for i in range(len(sdssobjects))])
    wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    flux = np.array([sdss_data_files[i][:,1] for i in range(len(sdssobjects))])
    fluxerr = np.array([sdss_data_files[i][:,2] for i in range(len(sdssobjects))])
    bp_map = np.array([sdss_data_files[i][:,3] for i in range(len(sdssobjects))])
    flux_cont = np.array([sdss_data_files[i][:,6] for i in range(len(sdssobjects))])



    # Generate SDSS composite
    wl_sdss = np.array([sdss_data_files[i][:,4] / (1 + redshifts[i]) for i in range(len(sdssobjects))])
    flux_sdss = np.array([sdss_data_files[i][:,5] for i in range(len(sdssobjects))])



    # Interpolate to a common wavelength:
    short_sdss = []
    tall_sdss = []
    for i in wl_sdss:
        short_sdss.append(min(i[np.where(i != 0)]))
        tall_sdss.append(max(i))
    short_sdss = min(short_sdss)
    tall_sdss = max(tall_sdss)

    step_sdss = 1.0 #CDELT
    wl_new_sdss = np.arange(short_sdss, tall_sdss, step_sdss)


    flux_new_sdss = np.zeros((len(redshifts),len(wl_new_sdss)))

    from gen_methods import smooth,medfilt

    # for n_sdss in range(np.shape(wl_sdss)[0]):
    #
    #
    #     flux_sdss[n_sdss] = (flux_sdss[n_sdss])[np.where(wl_sdss[n_sdss] != 0)]
    #     wl_sdss[n_sdss] = (wl_sdss[n_sdss])[np.where(wl_sdss[n_sdss] != 0)]
    #
    #
    #     # pl.plot(wl_sdss[n_sdss], flux_sdss[n_sdss], lw=0.2)
    #     # pl.semilogy()
    #     # pl.show()
    #
    #
    #
    #     #de-reddening
    #     flux_sdss[n_sdss] = cardelli_reddening(wl_sdss[n_sdss], flux_sdss[n_sdss], ebv[n_sdss])
    #
    #     #Interpolate
    #     f = interpolate.interp1d(wl_sdss[n_sdss],flux_sdss[n_sdss],kind='linear',bounds_error = False, fill_value=0.)
    #     # mask = (wl_new_sdss > 3020) & (wl_new_sdss < 3100)
    #     mask = (wl_new_sdss > 2400) & (wl_new_sdss < 2500)
    #     print(np.median(f(wl_new_sdss)[mask]))
    #     print(f(wl_new_sdss))
    #     norm_sdss = np.median(f(wl_new_sdss)[mask])
    #     flux_new_sdss[n_sdss] = f(wl_new_sdss)/norm_sdss
    #     # pl.plot(wl_new_sdss, flux_new_sdss[n_sdss], lw=0.2)
    #     # pl.semilogy()
    #     # pl.show()
    wl_sdss = wl_new_sdss
    flux_sdss = flux_new_sdss

    mean_sdss = np.zeros(np.shape(wl_sdss))
    for i_sdss, k_sdss in enumerate(flux_sdss.transpose()):

        #Mean

        # print(k_sdss)
        mean_sdss[i_sdss] = np.mean(k_sdss)
        # print(np.mean(k_sdss))











    # TODO Implement this interpolation as a more general method
    # Interpolate to a common wavelength:
    short = []
    tall = []
    for i in wl:
        short.append(min(i))
        tall.append(max(i))
    short = min(short)
    tall = max(tall)

    step = 0.2 #CDELT
    wl_new = np.arange(short, tall, step)
    print(short, tall, len(wl_new))
    flux_new = np.zeros((len(redshifts),len(wl_new)))
    flux_cont_new = np.zeros((len(redshifts),len(wl_new)))

    fluxerr_new = np.zeros((len(redshifts),len(wl_new)))
    bp_map_new = np.zeros((len(redshifts),len(wl_new)))
    from gen_methods import smooth,medfilt
    for n in range(np.shape(wl)[0]):
        #de-reddening
        flux[n] = cardelli_reddening(wl_obs[n], flux[n], ebv[n])
        flux_cont[n] = cardelli_reddening(wl_obs[n], flux_cont[n], ebv[n])

        fluxerr[n] = cardelli_reddening(wl_obs[n], fluxerr[n], ebv[n])
        #Interpolate
        f = interpolate.interp1d(wl[n],flux[n],kind='linear',bounds_error = False, fill_value=0.)
        j = interpolate.interp1d(wl[n],flux_cont[n],kind='linear',bounds_error = False, fill_value=0.)
        g = interpolate.interp1d(wl[n],fluxerr[n],kind='linear',bounds_error = False, fill_value=1.)
        h = interpolate.interp1d(wl[n],bp_map[n],kind='linear',bounds_error = False, fill_value=1.)
        mask = (wl_new > 7000) & (wl_new < 7500)
        # mask = (wl_new > 6900) & (wl_new < 7000)
        # mask = (wl_new > 2400) & (wl_new < 2500)
        # mask = (wl_new > 3020) & (wl_new < 3100)
        norm = np.median(f(wl_new)[mask])
        flux_new[n] = f(wl_new)/norm
        flux_cont_new[n] = j(wl_new)/norm
        fluxerr_new[n] = g(wl_new)/norm
        bp_map_new[n] = h(wl_new)
    #     pl.plot(wl_new, medfilt(flux_new[n], 1), lw=0.2)
    #     pl.plot(wl_new, medfilt(flux_cont_new[n], 1), lw=0.2)
    #
    # pl.semilogy()
    # pl.show()
    wl = wl_new
    flux = flux_new
    flux_cont = flux_cont_new
    fluxerr = fluxerr_new
    bp_map = bp_map_new


    for i, k in enumerate(bp_map):
        mask_cont = (wl < 1216)

        (bp_map[i])[mask_cont] = 0







    #------------------------- Combination -------------------------
    # TODO Methodize this to avoid microchangning
    # Weighted average:
    wmean = np.zeros(np.shape(wl))
    wmean_cont = np.zeros(np.shape(wl))
    mean = np.zeros(np.shape(wl))
    geo_mean = np.zeros(np.shape(wl))
    median = np.zeros(np.shape(wl))
    errofmean = np.zeros(np.shape(wl))
    errofwmean = np.zeros(np.shape(wl))
    std = np.zeros(np.shape(wl))
    import scipy
    import scikits.bootstrap as bootstrap
    CI_low = np.zeros(np.shape(mean))
    CI_high = np.zeros(np.shape(mean))
    for i, k in enumerate(flux.transpose()):
        mask = np.where(bp_map.transpose()[i] == 0)
        if len(k[mask]) != 0:
            #Weighted mean
            weight = 1./(np.array(fluxerr.transpose()[i][mask])**2)
            wmean[i] = np.average(k[mask], axis = 0, weights = weight)
            wmean_cont[i] = np.average((flux_cont.transpose()[i])[mask], axis = 0, weights = weight)
            errofwmean[i] = np.sqrt(1./np.sum(np.array(fluxerr.transpose()[i][mask])**-2.,axis=0))


            #Mean
            mean[i] = np.mean(k[mask])
            errofmean[i] = np.sqrt(np.sum((fluxerr.transpose()[i][mask])**2))


            #Geometric mean
            from scipy.stats.mstats import gmean
            geo_mean[i] = gmean(k[mask])

            #Median
            median[i] = np.median(k[mask])
            # Bootstrapping to get confidence intervals
            # if len(k[mask]) > 1:
            #     CI_low[i], CI_high[i] = bootstrap.ci(data=k[mask], statfunction=np.median, n_samples=100, method='bca')
            # else:
            #     CI_low[i], CI_high[i] = np.median(k[mask]), np.median(k[mask])
            # if (i % 1000==0):
            #     print(i)


            std[i] = np.std(flux.transpose()[i][mask])
        else:
            mean[i] = 0
            errofmean[i] = 0
            wmean[i] = 0
            wmean_cont[i] = 0
            errofwmean[i] = 0
            geo_mean[i] = 0
            median[i] = 0
            CI_low[i], CI_high[i] = 0, 0
            std[i] = 0


    #Insert zeros
    wl_sdss = np.concatenate([wl_sdss,np.zeros(len(wl) - len(wl_sdss))])
    mean_sdss = np.concatenate([mean_sdss,np.zeros(len(mean) - len(mean_sdss))])









    # TODO Methodize this?
    # Calculating the number of spectra that goes into the composite and appending to n_spec.
    n_spec = np.zeros(np.shape(mean))
    spec = []
    for i, n in enumerate(flux.transpose()):
        mask = np.where(bp_map.transpose()[i] == 0)
        n_spec[i] = len((n[mask])[np.where((n[mask]) != 0)])
        if len(n[np.where(n != 0)]) == len(redshifts):
            spec.append(n / np.median(n))

    std = std / np.sqrt(n_spec)




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
          ("std", np.float64), ("wl_sdss", np.float64), ("mean_sdss", np.float64), ("wmean_cont", np.float64) ]
    data = np.array(zip(wl, mean, errofmean, wmean, errofwmean, geo_mean, median, n_spec, std, wl_sdss, mean_sdss, wmean_cont ), dtype=dt)
    file_name = "Composite"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl mean mean_error wmean wmean_error geo_mean median n_spec std wl_sdss mean_sdss wmean_cont")#, fmt = ['%5.1f', '%2.15E'] )

    #Saving to .dat file
    dt = [("wl", np.float64), ("wmean", np.float64) ]
    data = np.array(zip(wl, wmean), dtype=dt)
    file_name = "XSH-Composite"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl wmean")#, fmt = ['%5.1f', '%2.15E'] )

    #Saving to .dat file
    dt = [("wl", np.float64), ("wmean_cont", np.float64) ]
    data = np.array(zip(wl, wmean_cont), dtype=dt)
    file_name = "XSH-Composite_cont"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl wmean_cont")#, fmt = ['%5.1f', '%2.15E'] )

if __name__ == '__main__':
    main()