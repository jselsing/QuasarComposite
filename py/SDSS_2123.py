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
    """Main script to prepare x-shooter observations for combination"""
    from matplotlib import rc_file
    rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    import plotting.plotting as plot
    from xshoo.combine import inter_arm_cut
    from scipy.interpolate import splrep,splev

    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/SDSS2123-0050'
    object_files = glob.glob(root_dir+'/OBJECT/*/*/*IDP*.fits')
    respose_files = glob.glob(root_dir+'/M.X*.fits')



    arms = ['UVB', 'VIS', 'NIR']


    for n in arms:
        print('In arm: '+n)


        #Read in object spectrum
        obser = [k for k in object_files if n in k]
        ob = fits.open(obser[0])
        wl = 10.0*ob[1].data.field('WAVE')[0]
        flux = ob[1].data.field('FLUX')[0]
        err = ob[1].data.field('ERR')[0]
        print(fits.open(k)[0].header['RA'])
        print(fits.open(k)[0].header['DEC'])

        #Read in master response curve
        master_response = [k for k in respose_files if n in fits.open(k)[0].header['HIERARCH ESO SEQ ARM']]
        resp = fits.open(master_response[0])
        response_wl = resp[1].data.field('LAMBDA')*10.0
        response = (resp[1].data.field('RESPONSE'))
        response_wl = response_wl
        interp = splrep(response_wl,response)

        #Apply master response function
        flux *= splev(wl,interp)




        pl.plot(wl,flux)
        pl.show()


http://skyserver.sdss.org/dr10/en/tools/explore/summary.aspx?id=0x112d06c121100010&spec=0x0a64188878006800&apid=




if __name__ == '__main__':
    main()