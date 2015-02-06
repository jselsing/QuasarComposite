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



if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    from scipy.interpolate import splrep,splev
    import plotting.plotting as plot
    from xshoo.combine import inter_arm_cut

    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/')
    sdssobjects = ['SDSS0021+0043/', 'SDSS0022+0124/', 'SDSS0820+1306/', 'SDSS1013+0851/', 'SDSS1101+0548/', 'SDSS1150-0023/', 'SDSS1219-0100/', 'SDSS1236-0331/', 'SDSS1249-0559/', 'SDSS1354-0013/', 'SDSS1431+0535/', 'SDSS1437-0147/', 'SDSS2123-0050/', 'SDSS2313+0034/']
    sdssobjects = ['SDSS0820+1306/', 'SDSS1150-0023/', 'SDSS1219-0100/', 'SDSS1236-0331/', 'SDSS1354-0013/', 'SDSS1431+0535/', 'SDSS1437-0147/']
    sdssobjects = ['SDSS0022+0124/', 'SDSS1013+0851/', 'SDSS1101+0548/', 'SDSS1249-0559/']
    # sdssobjects = ['SDSS1249-0559/']
    # arms = ['UVB', 'VIS', 'NIR']
    arms = ['VIS', 'NIR']
    for x,i in enumerate(sdssobjects):
        obs = glob.glob(root_dir+i+'*TELL*/*/*IDP*')
        print(obs)
        master_response = glob.glob(root_dir+i+'M.X*.fits')

        obj_name = i[-14:-1]
        wl_out = []
        flux_out = []
        test = []
        err_out = []
        start = []
        end = []
        print(obj_name)
        #arms = ['NIR']
        comb = []
        for l,n in enumerate(arms):
            print('In arm: '+n)
            obser = [k for k in obs if n in k]
            ob = fits.open(obser[0])
            wl = 10.0*ob[1].data.field('WAVE')[0]
            flux = ob[1].data.field('FLUX')[0]
            err = ob[1].data.field('ERR')[0]



            resp = fits.open(master_response[l+1])
            response_wl = resp[1].data.field('LAMBDA')*10.0
            # mask = np.where(response_wl < max(wl))[0]
            response = (resp[1].data.field('RESPONSE'))
            response_wl = response_wl
            interp = splrep(response_wl,response)
            response_interp = splev(wl,interp)






            print(len(response))
            mask = np.where(response_wl < max(wl))[0]
            print(len(response_wl[mask]))




            # pl.plot(response_wl, response)
            # pl.plot(wl, response_interp)
            pl.plot(wl,flux*response_interp, lw = 0.5, linestyle='steps-mid')
            pl.show()

