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

    wl = np.array([(sdss_data_files[i][:,0] / (1 + redshifts[i])) for i in range(len(sdssobjects))])
    wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    flux = np.array([(sdss_data_files[i][:,1] * (1 + redshifts[i])) for i in range(len(sdssobjects))])
    fluxerr = np.array([(sdss_data_files[i][:,2] * (1 + redshifts[i])) for i in range(len(sdssobjects))])
    bp_map = np.array([sdss_data_files[i][:,3] for i in range(len(sdssobjects))])
    flux_cont = np.array([sdss_data_files[i][:,6] * (1 + redshifts[i]) for i in range(len(sdssobjects))])
    n_obj = len(obj_list)






    filter = glob.glob('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SLOAN_SDSS.i.dat')[0]
    filter = np.genfromtxt(filter)
    wl_filt = filter[:,0]
    filt = filter[:,1]


    flux_johan = []
    # fig, ax = pl.subplots(2, sharex=True)
    # for mask_ran in range(1500,9000, 1500):
    pow_slope = []
    indi_pow = []
    # norm_reg = [1420, 3000, 4000, 5600, 6800]
    # norm_reg = np.arange(1420, 8000, 500)
    norm_reg = [6800]
    # print(norm_reg)
    for mask_ran in norm_reg:
        flux_0 = flux.copy()
        flux_cont_0 = flux_cont.copy()
        fluxerr_0 = fluxerr.copy()


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
        filt_new = np.zeros((n_obj,n_wl))

        for n in range(n_obj):
            #de-reddening
            flux_0[n] = cardelli_reddening(wl_obs[n], flux_0[n], ebv[n])
            flux_cont_0[n] = cardelli_reddening(wl_obs[n], flux_cont_0[n], ebv[n])
            fluxerr_0[n] = cardelli_reddening(wl_obs[n], fluxerr_0[n], ebv[n])

        for n in range(n_obj):
            #Interpolate
            flux_new[n] = common_wavelength(wl[n], wl_new, flux_0[n])
            flux_cont_new[n] = common_wavelength(wl[n], wl_new, flux_cont_0[n])
            fluxerr_new[n] = common_wavelength(wl[n], wl_new, fluxerr_0[n], fill_value=1.0)
            bp_map_new[n] = common_wavelength(wl[n], wl_new, bp_map[n], fill_value=1.0)


        # print(np.shape(flux))
        flux_cont_new = np.vstack((flux_cont_new , 1e-15 * (wl_new/10000.0)**(-1.7)))
        redshifts = np.concatenate((redshifts, [1.2]))
        obj_list.append('Pure power law')


        from astropy.cosmology import FlatLambdaCDM
        from astropy.cosmology import Planck13 as cosmo
        # cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        filt_new =  common_wavelength(wl_filt, wl_new, filt, fill_value=0.0)
        miz0 = []
        for n in range(n_obj+1):
            fig, ax = pl.subplots(1)
            ax.plot(wl_new, filt_new)
            ax2 = ax.twinx()
            ax2.plot(wl_new, medfilt(flux_cont_new[n],31))
            pl.show()
            prod = medfilt(filt_new * flux_cont_new[n], 29)


            numerator = np.sum(prod * wl_new)
            denom = np.sum(filt_new * (3e18/wl_new))
            f_nu = numerator / denom
            i_band_mag = -2.5 * np.log10(f_nu) - 48.6
            dl = (cosmo.luminosity_distance(redshifts[n])) * 1e5
            M = -5 * np.log10(dl.value) + i_band_mag
            miz0.append(M)
            print(obj_list[n])
            print(M, i_band_mag, redshifts[n])
        # pl.plot(wl_new, prod)
        # pl.show()

        # print('')

        wl_new_shift = wl_new * (1+2)
        filt_new =  common_wavelength(wl_filt, wl_new_shift, filt, fill_value=0.0)
        miz2 = []
        for n in range(n_obj+1):
            flux_cont_new[n] /= (1+2)
            prod = medfilt(filt_new * flux_cont_new[n], 29)

            numerator = np.sum(prod * wl_new_shift)
            denom = np.sum(filt_new * (3e18/wl_new_shift))
            f_nu = numerator / denom
            i_band_mag = -2.5 * np.log10(f_nu) - 48.6
            dl = (cosmo.luminosity_distance(redshifts[n])) * 1e5
            # dl = 1
            M = -5 * np.log10(dl.value) + i_band_mag
            miz2.append(M)
            print(obj_list[n])
            print(M, i_band_mag, redshifts[n])

        print(np.array(miz2) - np.array(miz0))

        # pl.plot(wl_new_shift, prod)
        # pl.show()



        for n in range(n_obj):
            #Normalise
            mask = (wl_new > mask_ran) & (wl_new < mask_ran + 100)
            norm = np.median(flux_new[n][mask])
            flux_new[n] /= norm
            flux_cont_new[n] /= norm
            fluxerr_new[n] /= norm

            # flux_cont_new[np.where(flux_cont_new == np.NaN)] = 1

            par_guess = [2e6, -1.5]
            # mask = (wl_new > 1350) & (wl_new < 1365) | (wl_new > 4200) & (wl_new < 4230) | (wl_new > 5500) & (wl_new < 6035) #| (wl_new > 7800) & (wl_new < 7950)
            mask = (wl_new > 1300) & (wl_new < 1350) | (wl_new > 1425) & (wl_new < 1475) | (wl_new > 5500) & (wl_new < 5800) | (wl_new > 7300) & (wl_new < 7500)
            popt, pcov = optimize.curve_fit(power_law, wl_new[mask], flux_cont_new[n][mask], p0=par_guess,
                                            sigma=fluxerr_new[n][mask] , absolute_sigma=True, maxfev=2000)
            # print(*popt)
            indi_pow.append(popt[1])


        #     ax.plot(wl_new, power_law(wl_new, *popt), '--', lw = 0.5)
        #     ax.plot(wl_new, medfilt(flux_cont_new[n], 5), lw = 0.5, label= str(mask_ran))
        # pl.semilogx()
        # pl.semilogy()
        # pl.show()

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
        for i, k in enumerate(flux_cont_new.transpose()):
            mask = np.where(bp_map_new.transpose()[i] == 0)
            # mask = np.where(bp_map_new.transpose()[i] > -10000)
            # hej = True
            # if hej:
            if len(k[mask]) != 0:
                #Weighted mean
                e = np.array(fluxerr_new.transpose()[i][mask])
                weight = 1. / e ** 2
                wmean[i] = np.average(flux_new.transpose()[i][mask], axis = 0, weights = weight)
                wmean_cont[i] = np.average(k[mask], axis = 0, weights = weight)
                errofwmean[i] = np.sum(weight,axis=0) ** -0.5

                #Mean
                mean[i] = np.mean(k[mask])
                errofmean[i] = np.sqrt(np.sum((fluxerr_new.transpose()[i][mask])**2))

                #Geometric mean
                from scipy.stats.mstats import gmean

                geo_mean[i] = gmean((k[mask]))


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


                std[i] = np.std(k[mask])
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
        # test =  wmean_cont / norm
        # ii  = (wl_new > 6900) & (wl_new < 6900 + 50)
        # norm = np.median(geo_mean[ii])
        # geo_mean_norm = geo_mean / norm
        # ax[0].plot(wl_new, medfilt(geo_mean_norm, 5), lw = 0.5, label= str(mask_ran) + '_geo')


        # ii  = (wl_new > 6900) & (wl_new < 6900 + 50)
        # norm = np.median(wmean_cont[ii])
        # wmean_cont_norm = wmean_cont / norm
        # ax[0].plot(wl_new, medfilt(wmean_cont / norm, 5), lw = 0.5, label= str(mask_ran) + '_wmean')

        # ii  = (wl_new > 6900) & (wl_new < 6900 + 50)
        # norm = n.median(median[ii])
        # ax[0].plot(wl_new, 1 - medfilt(wmean_cont_norm, 5)/medfilt(median/norm, 5), lw = 0.5, label= str(mask_ran) + '_median')
        # test = 1 - medfilt(wmean_cont_norm[median != 0.],5 ) /medfilt((median[median != 0.]/norm),5)
        # print(test)
        # print("""{0}  +- {1} """.format(np.mean(test), np.std(test)))

        # ii  = (wl_new > 6900) & (wl_new < 6900 + 50)
        # norm = np.median(mean[ii])
        # ax[0].plot(wl_new, 1 - medfilt(wmean_cont_norm, 5)/medfilt(mean / norm, 5), lw = 0.5, label= str(mask_ran) + '_mean')
        # test = 1 - medfilt(wmean_cont_norm[mean != 0.],5 )/medfilt((mean[mean != 0.] / norm),5 )
        # print("""{0}  +- {1} """.format(np.mean(test), np.std(test)))
        par_guess = [1, -1.0]
        par_guess2 = [1, 3200, -1.0, -1.46]
        # par_guess3 = [1, -1.0, -0.00001]



        # pl.plot(wl, power_law2(wl, par_guess))
        # pl.show()
        # geo_mean_norm[np.where(np.isnan(geo_mean_norm) == True)] = 0

        # mask = (wl_new > 1350) & (wl_new < 1365) | (wl_new > 4200) & (wl_new < 4230) | (wl_new > 5500) & (wl_new < 6035) | (wl_new > 7800) & (wl_new < 7950)
        mask = (wl_new > 1300) & (wl_new < 1350) | (wl_new > 1425) & (wl_new < 1475) | (wl_new > 5500) & (wl_new < 5800) | (wl_new > 7300) & (wl_new < 7500)
        popt_geo, pcov_geo = optimize.curve_fit(power_law, wl_new[mask], geo_mean[mask], p0=par_guess)
        # popt_wmean, pcov_wmean = optimize.curve_fit(power_law, wl_new[mask], wmean_cont[mask], p0=par_guess,
        #                                             sigma=errofwmean[mask], absolute_sigma=True, maxfev = 2000)
        popt_wmean, pcov_wmean = optimize.curve_fit(power_law, wl_new[mask], wmean_cont[mask], p0=par_guess, maxfev = 2000)
        # popt_mean, pcov_mean = optimize.curve_fit(power_law, wl_new[mask], mean[mask], p0=par_guess,
        #                                           sigma=errofmean[mask], absolute_sigma=True, maxfev = 2000)
        popt_mean, pcov_mean = optimize.curve_fit(power_law, wl_new[mask], mean[mask], p0=par_guess, maxfev = 2000)
        popt_median, pcov_median = optimize.curve_fit(power_law, wl_new[mask], mean[mask], p0=par_guess)
        print("""Composite fit slope geo...{0} +- {1}""".format(popt_geo[1], np.diag(pcov_geo)[1]))
        print("""Composite fit slope wmean...{0} +- {1}""".format(popt_wmean[1], np.diag(pcov_wmean)[1]))
        print("""Composite fit slope mean...{0} +- {1}""".format(popt_mean[1], np.diag(pcov_mean)[1]))
        print("""Composite fit slope median...{0} +- {1}""".format(popt_median[1], np.diag(pcov_median)[1]))
        print("""Individual slope mean...{0} +- {1}""".format(np.mean(indi_pow), np.std(indi_pow)))
        print("""Individual slope median...{0} +- {1}""".format(np.median(indi_pow), np.std(indi_pow)))

        # ax.plot(wl_new, power_law2(wl_new, *popt2) , '--', lw = 0.5)
        # ax[0].plot(wl_new, power_law(wl_new, *popt_wmean) , 'b--', lw = 0.5)
        # ax[0].plot(wl_new, power_law(wl_new, *popt_mean) , 'b--', lw = 0.5)

        pow_slope.append(popt[1])

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



        # flux_johan.append(geo_mean.transpose())

    # np.savetxt('test2.dat', zip(wl_new, *flux_johan), header="wl 1500 3000 4500 6000 7500")#, fmt = ['%5.1f', '%2.15E'] )

    # ax[0].semilogy()
    # ax[0].semilogx()
    # ax[1].semilogy()
    # ax[1].semilogx()
    # ax[0].legend()
    # pl.show()



    # fig, ax = pl.subplots(1)
    # ax.plot(norm_reg, pow_slope, '.')
    # print("""Power law slopes:
    #     Mean = {0}
    #     Std  = {1}
    # """.format(np.mean(pow_slope), np.std(pow_slope)))
    # ax.set_xlabel('Normalisation region')
    # ax.set_ylabel('Power law slope')
    # ax.invert_yaxis()
    # pl.show()

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