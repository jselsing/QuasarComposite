

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
import matplotlib.pylab as pl
import seaborn as sns; sns.set_style('ticks')
#cmap = sns.cubehelix_palette(n_colors=6, start=1, rot=0.2, gamma=1.0, hue=0.8, light=0.85, dark=0.15, reverse=True, as_cmap=False)
cmap = sns.color_palette("cubehelix", 6)

import matplotlib

from math import sqrt
SPINE_COLOR = 'gray'

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES
    size = 16
    params = {'backend': 'Qt4Agg',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': size, # fontsize for x and y labels (was 10)
              'axes.titlesize': size,
              'font.size': size, # was 10
              'legend.fontsize': size, # was 10
              'xtick.labelsize': size,
              'ytick.labelsize': size,
              #'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)


def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(True)

    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax

def main():
    """Main script to prepare x-shooter observations for combination"""
    # from matplotlib import rc_file
    # rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
    from astropy.io import fits
    import glob
    # import matplotlib.pyplot as pl
    latexify()

    import numpy as np
    # import plotting.plotting as plot
    from xshoo.combine import inter_arm_cut
    from scipy.interpolate import splrep,splev

    #Files
    obj_name = 'SDSS0820+1306'
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'+obj_name
    object_files = glob.glob(root_dir+'/OBJECT/*IDP*.fits')
    # respose_files = glob.glob(root_dir+'/M.X*.fits')
    transmission_files = glob.glob(root_dir+'/transmission*.fits')
    arms = ['UVB', 'VIS', 'NIR']
    wl_out = []
    flux_out = []
    flux_uncorr_out = []
    err_out = []
    start = []
    end = []

    for n in arms:
        print('In arm: '+n)


        #Read in object spectrum
        obser = [k for k in object_files if n in k]
        ob = fits.open(obser[0])
        wl = 10.0*ob[1].data.field('WAVE')[0]
        flux = ob[1].data.field('FLUX')[0]
        err = ob[1].data.field('ERR')[0]

        print(np.shape(flux))
        #Read in master response curve
        # master_response = [k for k in respose_files if n in fits.open(k)[0].header['HIERARCH ESO SEQ ARM']]
        # resp = fits.open(master_response[0])
        # response_wl = resp[1].data.field('LAMBDA')*10.0
        # response = (resp[1].data.field('RESPONSE'))
        # response_wl = response_wl
        # interp = splrep(response_wl,response)

        #Apply master response function
        # flux *= splev(wl,interp)
        # err *= splev(wl,interp)

        wl_tmp, flux_uncorr, err_tmp, start_tmp, end_tmp = inter_arm_cut(wl, flux, err, n, start, end)
        if n== 'VIS' or n== 'NIR':
            transmission = fits.open([k for k in transmission_files if n in k][0])[0].data
            print(np.shape(transmission))
            for j, k in enumerate(transmission):
                if k <= 1e-10:
                    transmission[j] = 1

            flux /= transmission
            err /= transmission

        wl, flux, err, start, end = inter_arm_cut(wl, flux, err, n, start, end)


        wl_out.append(wl)
        flux_out.append(flux)
        err_out.append(err)
        flux_uncorr_out.append(flux_uncorr)


    wl_out = np.hstack(wl_out)
    flux_out = np.hstack(flux_out)
    err_out = np.hstack(err_out)
    flux_uncorr_out = np.hstack(flux_uncorr_out)


    bp_map = []
    npix = 0
    for j , (k, l) in enumerate(zip(flux_out[:-1],err_out[:-1])):
        if k > 1.1 * flux_out[j-1] or k < 0:
           bp_map.append(1)
           npix += 1
        elif k < 0.90 * flux_out[j-1] or k < 0:
           bp_map.append(1)
           npix += 1
        else:
           bp_map.append(0)
    bp_map.append(1)
    print(npix)
    print(len(bp_map))
    print((len(bp_map) - npix/len(bp_map))/len(bp_map))


    import json
    import urllib2

    query_terms = dict()
    query_terms["ra"] = str(ob[0].header['RA'])+'d' #"185.1d"
    query_terms["dec"] = str(ob[0].header['DEC'])  #"56.78"
    query_terms["radius"] = "5.0"


    url = "http://api.sdss3.org/spectrumQuery?" + '&'.join(["{0}={1}".format(key, value) for key, value in query_terms.items()])
    print(url)
    # make call to API
    response = urllib2.urlopen(url)

    # read response, converting JSON format to Python list
    matching_ids = json.loads(response.read())
    print(json.dumps(matching_ids, indent=4))

    # get the first id
    spec_id = matching_ids[0]

    url = "http://api.sdss3.org/spectrum?id={0}&format=json".format(spec_id)

    response = urllib2.urlopen(url)
    result = json.loads(response.read())
    SDSS_spectrum = result[spec_id]

    wl_sdss = np.array(SDSS_spectrum["wavelengths"])
    flux_sdss =  np.array(SDSS_spectrum["flux"])
    z_sdss = (np.array(SDSS_spectrum["z"]))
    z_sdss_err = ((np.array(SDSS_spectrum["z_err"])))


    wl_sdss = np.concatenate([wl_sdss,np.zeros(len(wl_out) - len(wl_sdss))])
    flux_sdss = np.concatenate([flux_sdss,np.zeros(len(flux_out) - len(flux_sdss))])


    # Load linelist
    fit_line_positions = np.genfromtxt('data/fitlinelist.txt', dtype=None)
    linelist = []
    for n in fit_line_positions:
        linelist.append(n[1])
    linelist = np.array(linelist)


    from methods import wavelength_conversion
    linelist = wavelength_conversion(linelist, conversion='vacuum_to_air')


    #Cut out fitting region
    mask = np.logical_and(wl_out > 10570, (wl_out < 10700))
    wl_fit = wl_out[mask]
    flux_fit = flux_out[mask]
    fluxerr_fit = err_out[mask]

    fluxerr_new = []
    for j, (k, l) in enumerate(zip(flux_fit,fluxerr_fit)):
        if k > 1.5 * flux_fit[j-2] and k > 0:
            fluxerr_new.append(l*50)
        elif k < 0.75 * flux_fit[j-2] and k > 0:
            fluxerr_new.append(l*50)
        else:
            fluxerr_new.append(l)
    from gen_methods import smooth
    fluxerr_fit = smooth(np.array(fluxerr_new), window_len=15, window='hanning')


    #Fit continuum and subtract
    from methods import continuum_fit
    from numpy.polynomial import chebyshev
    cont, chebfit = continuum_fit(wl_fit, flux_fit, fluxerr_fit, edge_mask_len=20)
    chebfitval = chebyshev.chebval(wl, chebfit)


    #Define models to use
    from methods import voigt,gauss
    def model1(t,  amp2, sig22g, sig22l, z):
            tmp = voigt(t, abs(amp2), (1+z)*linelist[2], sig22g, sig22l)
            return tmp

    def model2(t, amp2, sig22g, z):
            tmp = gauss(t, abs(amp2), (1+z)*linelist[2], sig22g)
            return tmp


    #Initial parameters
    init_vals = [6e-17,10, z_sdss]
    y_fit_guess = model2(wl_fit, *init_vals) + cont


    #Fit
    import scipy.optimize as op
    np.random.seed(12345)
    y_op = []
    vals = []
    for i in np.arange(10000):
        print('Iteration: ', i)
        resampled_spec = np.random.normal(flux_fit, fluxerr_fit)

        cont, chebfit = continuum_fit(wl_fit, resampled_spec, fluxerr_fit, edge_mask_len=20)
        chebfitval = chebyshev.chebval(wl, chebfit)

        best_vals, covar = op.curve_fit(model2, wl_fit, resampled_spec - cont, sigma=fluxerr_fit, absolute_sigma=True, p0=init_vals)
        vals.append(best_vals)


    up = (np.percentile(vals, 84, axis = 0)[2] - np.mean(vals, axis = 0)[2])
    down = (np.percentile(vals, 16, axis = 0)[2] - np.mean(vals, axis = 0)[2])



    v_bary = ob[0].header['HIERARCH ESO QC VRAD BARYCOR']
    c_km = (2.99792458e8/1000.0)

    print("""Curve_fit results:
        Redshift = {0} + {1} - {2} (SDSS: {3} +- {4})
    """.format(np.mean(vals, axis = 0)[2] + v_bary /c_km, up, down, z_sdss, z_sdss_err))


    # #Calculate best fit values + confidence values
    # y_op = model2(wl_fit, *best_vals) + cont
    # best_vals[2] = z_sdss
    # y_sdss = model2(wl_fit, *best_vals) + cont
    #
    # from methods import ConfInt
    # y_op_lower, y_op_upper = ConfInt(wl_fit, model2, best_vals, covar, [16,84]) + cont



    # #Overplot lines
    # fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)
    # linelist = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    # linelist = np.array(linelist)

    # # print((1+z_op)*linelist[1])
    # mask = (wl_out < (1+z_op)*linelist[1])
    # wave = wl_out[mask]
    # flux = flux_out[mask]
    # import continuum_mark.interactive
    # normalise = continuum_mark.interactive.continuum_mark(wl_out[mask], flux_out[mask], err_out[mask])
    # normalise.endpoint = 'n' #str(raw_input('Insert endpoint before interpolation(y/n)? '))
    #
    # normalise.run()
    # pl.show()
    # cont_out = np.concatenate([normalise.continuum,flux_out[~mask]])




    # # Plotting for the paper
    # from plotting import plotting

    # fig = plotting.plot_data(wl_out, flux_out, xrng=xrng, yrng=yrng, lw=0.8,
    #         label='SDSS0820+1306', label_error='Error',
    #         plot_top_tick=True, z=z_op)
    # ax1, ax2, ax3 = fig.axes


    # wl_fit /= (1+z_op)
    # wl_out /= (1+z_op)
    # y_op *= (1+z_op)
    # y_op_lower *= (1+z_op)
    # y_op_upper *= (1+z_op)
    # flux_out *= (1+z_op)
    #
    #
    #
    # ratio = (1.0 + np.sqrt(5.0))/2.0
    # fig, ax = pl.subplots(figsize=(5*ratio, 5))
    # # fig, ax = pl.subplots()
    #
    #
    # # #Overplot lines
    # # for p in range(len(fit_line_positions)):
    # #     xcoord = linelist[p]*(1+z_op)
    # #     mask = (wl > xcoord - 1) & (wl < xcoord + 1)
    # #     y_val = np.mean(flux[mask])
    # #     pl.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.75)
    # #     pl.annotate(fit_line_positions[p,][0],xy=(xcoord, y_val * 1.15 ),fontsize='small')
    # #
    # #
    # ax.plot(wl_out, flux_out , linestyle = 'steps-mid', lw=0.6)
    # ax.plot(wl_out, err_out)
    # #pl.plot(wl_fit,flux_fit, color = 'green', lw = 1.0, alpha = 0.5)
    # #pl.plot(wl_fit, y_fit_guess)
    # ax.plot(wl_fit, y_op)
    # # ax.plot(wl_fit, y_sdss)
    # ax.plot(wl_fit, fluxerr_fit)
    #
    # ax.fill_between(wl_fit, y_op_lower, y_op_upper, alpha = 0.2)
    #
    # #Overplot lines
    # fit_line_positions = np.genfromtxt('data/plotlinelist.txt', dtype=None)
    # import lineid_plot
    #
    # linelist = []
    # linenames = []
    # for n in fit_line_positions:
    #     linelist.append(n[1])
    #     linenames.append(n[0])
    # print(linenames)
    # val = []
    # for p in range(len(linelist)):
    #     xcoord = linelist[p]
    #     mask = (wl_out > xcoord - 1) & (wl_out < xcoord + 1)
    #     try:
    #         y_val = np.mean(flux_out[mask])
    #         val.append(2 * y_val)
    #     except:
    #         pass
    # # arrow_tips = val, arrow_tip=arrow_tips
    # lineid_plot.plot_line_ids(wl_out, flux_out, linelist, linenames, ax=ax)
    # for i in ax.lines:
    #     if '$' in i.get_label():
    #         i.set_alpha(0.3)
    #
    #






    # pl.xlim((10100/ (1+z_op), 10800/ (1+z_op)))
    # pl.ylim((3.75e-16* (1+z_op), 7.4e-16* (1+z_op)))
    # ax.set_xlabel(r'Restframe Wavelength [$\AA$]')
    # ax.set_ylabel(r'Flux [erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$]')
    # # for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
    # #          ax.get_xticklabels(which='both') + ax.get_yticklabels(which='both')):
    # #     item.set_fontsize(16)
    #
    # format_axes(ax)
    #
    # pl.tight_layout()
    #
    # pl.savefig("../documents/figs/LineFit.pdf", clobber=True)
    # pl.show()
    flag = 1


    from astroquery.irsa_dust import IrsaDust
    import astropy.coordinates as coord
    import astropy.units as u
    C = coord.SkyCoord(ob[0].header['RA']*u.deg, ob[0].header['DEC']*u.deg, frame='fk5')
    dust_image = IrsaDust.get_images(C, radius=2 *u.deg, image_type='ebv')[0]
    ebv = np.mean(dust_image[0].data)
    print(ebv)




    #Saving to .dat file
    # dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64), ("bp map", np.float64),
    #       ("wl_sdss", np.float64), ("flux_sdss", np.float64) , ("flux_cont", np.float64) ]
    # data = np.array(zip(wl_out, flux_uncorr_out, err_out, bp_map, wl_sdss, flux_sdss, cont_out), dtype=dt)
    # file_name = "Telluric_uncorrected_science"
    # np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl flux fluxerror bp_map wl_sdss flux_sdss cont")#, fmt = ['%5.1f', '%2.15E'] )
    #
    #
    #
    # #Saving to .dat file
    # dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64), ("bp map", np.float64),
    #       ("wl_sdss", np.float64), ("flux_sdss", np.float64) , ("flux_cont", np.float64) ]
    # data = np.array(zip(wl_out, flux_out, err_out, bp_map, wl_sdss, flux_sdss, cont_out), dtype=dt)
    # file_name = "Telluric_corrected_science"
    # np.savetxt(root_dir+"/"+file_name+".dat", data, header="wl flux fluxerror bp_map wl_sdss flux_sdss cont")#, fmt = ['%5.1f', '%2.15E'] )
    #
    #
    # #Saving to .dat file
    # dt = [("z_op", np.float64), ("z_sdss", np.float64), ("flag", np.float64), ("ebv", np.float64)]
    # data = np.array(zip([z_op], [z_sdss], [flag], [ebv]), dtype=dt)
    # file_name = "Object_info"
    # np.savetxt(root_dir+"/"+file_name+".dat", data, header="z_op z_sdss flag ebv ") #, fmt = ['%5.1f', '%2.15E'] )





if __name__ == '__main__':
    main()
