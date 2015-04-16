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
    # redshifts = [1.1250, 1.9798, 1.5826, 1.8458, 1.5123, 2.0998, 1.3092]
    # ebv = [0.0253, 0.0201, 0.0279, 0.0268, 0.0330, 0.0274, 0.0378]
    # TODO Is there a better way to import the txt files?
    sdss_data_files = np.array([np.genfromtxt(i) for i in sdssobjects])

    # wl = np.array([sdss_data_files[i][:,0] / (1 + redshifts[i]) for i in range(len(sdssobjects))])
    # wl_obs = np.array([sdss_data_files[i][:,0] for i in range(len(sdssobjects))])
    # flux = np.array([sdss_data_files[i][:,1] for i in range(len(sdssobjects))])
    # fluxerr = np.array([sdss_data_files[i][:,2] for i in range(len(sdssobjects))])
    # bp_map = np.array([sdss_data_files[i][:,3] for i in range(len(sdssobjects))])
    # flux_cont = np.array([sdss_data_files[i][:,6] for i in range(len(sdssobjects))])



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

    for n_sdss in range(np.shape(wl_sdss)[0]):


        flux_sdss[n_sdss] = (flux_sdss[n_sdss])#[np.where(wl_sdss[n_sdss] != 0)]
        wl_sdss[n_sdss] = (wl_sdss[n_sdss])#[np.where(wl_sdss[n_sdss] != 0)]


        # pl.plot(wl_sdss[n_sdss], flux_sdss[n_sdss], lw=0.2)
        # pl.semilogy()
        # pl.show()



        #de-reddening
        flux_sdss[n_sdss] = cardelli_reddening(wl_sdss[n_sdss], flux_sdss[n_sdss], ebv[n_sdss])

        #Interpolate
        f = interpolate.interp1d(wl_sdss[n_sdss],flux_sdss[n_sdss],kind='linear',bounds_error = False, fill_value=0.)
        # mask = (wl_new_sdss > 3020) & (wl_new_sdss < 3100)
        mask = (wl_new_sdss > 2400) & (wl_new_sdss < 2500)
        norm_sdss = np.median(f(wl_new_sdss)[mask])
        flux_new_sdss[n_sdss] = f(wl_new_sdss)/norm_sdss
        # pl.plot(wl_new_sdss, flux_new_sdss[n_sdss], lw=0.2)
        # pl.semilogy()
        # pl.show()
    wl_sdss = wl_new_sdss
    flux_sdss = flux_new_sdss

    mean_sdss = np.zeros(np.shape(wl_sdss))
    for i_sdss, k_sdss in enumerate(flux_sdss.transpose()):

        #Mean
        mean_sdss[i_sdss] = np.mean(k_sdss)

    pl.plot(wl_sdss, mean_sdss)
    pl.show()

    #Saving to .dat file
    dt = [("wl_sdss", np.float64), ("wmean_sdss", np.float64) ]
    data = np.array(zip(wl_sdss, mean_sdss), dtype=dt)
    file_name = "SDSS-Composite"
    np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl_sdss wmean_sdss")#, fmt = ['%5.1f', '%2.15E'] )



if __name__ == '__main__':
    main()