#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


import seaborn as sns
import matplotlib.pyplot as pl
import numpy as np

def get_sdss_spectra(outfile = "outfile", N_spec = 5):
    from urllib2 import HTTPError
    from astroquery.sdss import SDSS
    # query = "SELECT TOP 1000 p.objid, p.dec, p.r,p.i, p.run, p.rerun, p.camcol, p.field, s.specobjid, s.class, s.z as redshift FROM PhotoObj AS p JOIN SpecObj AS s ON s.bestobjid = p.objid WHERE p.r BETWEEN 0 AND 17.0 AND s.class = 'QSO' AND s.z BETWEEN 1.0 AND 2.3 AND p.dec >= 15.0"
    query = "SELECT TOP "+str(N_spec)+" specObjID, plate, mjd, subClass, fiberID FROM SpecPhoto WHERE (class = 'QSO') AND" \
                                      " (psfmag_r <= 17.0) AND (dec >= 15.0) AND (z BETWEEN 1.0 AND 2.3) AND zwarning = 0 AND" \
                                      " (subClass = 'BROADLINE') AND nChild = 0 AND (mode = 1) AND ((0x10000000) != 0)" \
                                      " AND (bossprimary= 0) AND programname = 'legacy'"




    res = SDSS.query_sql(query)
# (subClass = 'BROADLINE') AND

    print(res['subClass'])
    spectra = []
    var = []
    waves = []
    mask = []
    z = []

    # print(res['plate'], res['mjd'], res['fiberID'])
    num_skipped = 0
    count = 1
    n_spec = len(res['specObjID'])
    for i in range(n_spec):
        # print(res['subClass'][i])
        try:
            sp = SDSS.get_spectra(plate=res['plate'][i], mjd=res['mjd'][i], fiberID=res['fiberID'][i])[0]
            data = (sp[1].data)

            wave = (10**data.field('loglam'))
            flux = data.field('flux')
            err = data.field('ivar')
            masking = data.field('and_mask')


            mask.append(masking)
            z.append(sp[2].data.field('Z'))
            spectra.append(flux)
            var.append(err)
            waves.append(wave)
            # print(res['plate'][i],res['mjd'][i], res['fiberID'][i])
            # pl.plot(wave, flux)
            # pl.show()
            count += 1
        except HTTPError:
            num_skipped += 1
            print("%i, %i, %i not found" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
            continue
        except ValueError:
            num_skipped += 1
            print("%i, %i, %i ValueError" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
            continue
        except TypeError:
            num_skipped += 1
            print("%i, %i, %i TypeError" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
            continue

        print('Number of spectrum processed: {0} out of {1}'.format(count, n_spec - num_skipped))
    print("   %i spectra skipped" % num_skipped)



    np.savez(outfile,
             wave = waves,
             spectra=spectra,
             var = var,
             mask = mask,
             plate = res['plate'],
             mjd = res['mjd'],
             fiberID = res['fiberID'],
             z = z
             )





def treat_sdss_spectra(outfile = "outfile"):

    data_in = np.load(outfile)
    spectra = data_in['spectra']
    var = data_in['var']
    mask = data_in['mask']
    waves = data_in['wave']
    redshifts = data_in['z']
    n_obj = len(redshifts)
    # print(redshifts)
    # print(np.shape(waves))
    # print(waves / (redshifts))
    for n, i in enumerate(redshifts):
        spectra[n] *= (1 + i)
        waves[n] /= (1 + i)

    short = []
    tall = []
    for wl_i in waves:
        short.append(min(wl_i))
        tall.append(max(wl_i))
    short = min(short)
    tall = max(tall)
    step = 2 #CDELT
    wl_new = np.arange(short, tall, step)
    n_wl = len(wl_new)
    spectra_new = np.zeros((n_obj,n_wl))
    var_new = np.zeros((n_obj,n_wl))

    from methods import common_wavelength
    for n in range(n_obj):
        #Interpolate
        spectra_new[n] = common_wavelength(waves[n], wl_new, spectra[n])
        var_new[n] = common_wavelength(waves[n], wl_new, var[n])

    norm_reg = 2850
    #fig, ax = pl.subplots()
    #Fitting power laws
    from scipy import optimize

    def power_law(x_tmp, a_tmp, k_tmp):

        tmp = a_tmp * x_tmp ** k_tmp
        return tmp






    indi_pow = []
    for n in range(n_obj):
        #Normalise
        mask = (wl_new > norm_reg) & (wl_new < norm_reg + 50)
        norm = np.median(spectra_new[n][mask])

        spectra_new[n] /= norm
        var_new[n] /= norm
        par_guess = [2e6, -1.5]
        mask = (wl_new > 1300) & (wl_new < 1350) | (wl_new > 1425) & (wl_new < 1475) | (wl_new > 5500) & (wl_new < 5800) | (wl_new > 7300) & (wl_new < 7500)
        popt, pcov = optimize.curve_fit(power_law, wl_new[mask], spectra_new[n][mask], p0=par_guess,
                                        sigma=var_new[n][mask] , absolute_sigma=True, maxfev=2000)

        indi_pow.append(popt[1])


    print("""Individual slope mean...{0} +- {1}""".format(np.mean(indi_pow), np.std(indi_pow)))
    print("""Individual slope median...{0} +- {1}""".format(np.median(indi_pow), np.std(indi_pow)))


     #   ax.plot(wl_new, spectra_new[n])
    #pl.show()









    wmean = np.zeros(np.shape(wl_new))
    errofwmean = np.zeros(np.shape(wl_new))
    for i, k in enumerate(spectra_new.transpose()):
        mask = np.where(k != 0.0)
        e = np.array(var_new.transpose()[i])[mask]
        weight = 1. / e ** 2
        wmean[i] = np.average(k[mask], axis = 0)#, weights = weight)
        errofwmean[i] = np.sum(weight,axis=0) ** -0.5


    pl.plot(wl_new, wmean)
    pl.loglog()
    pl.show()

    #Saving to .dat file
    dt = [("sdss_wl", np.float64), ("sdss_compo", np.float64)]
    data = np.array(zip(wl_new, wmean), dtype=dt)
    file_name = "sdss_compo"
    np.savetxt('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'+file_name+'.dat', data, header="wl flux") #, fmt = ['%5.1f', '%2.15E'] )


    # for i, k, l in zip(data_in['plate'], data_in['mjd'], data_in['fiberID']):
    #     print(str(i)+'+'+str(k)+'+'+str(l))



if __name__ == '__main__':

    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    file_name = "SDSS_spectra"
    outfile = root_dir+"/"+file_name+'.npz'

    #get_sdss_spectra(outfile = outfile, N_spec= 150)
    treat_sdss_spectra(outfile = outfile)

