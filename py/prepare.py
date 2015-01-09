# -*- coding: utf-8 -*-
"""
Call as Analyser.analysis(ebv,z,sdss,name)
"""
from __future__ import division, print_function

__all__ = []
__author__ = 'jselsing'

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')



# def cut(array_to_cut, boundary_array, lover_bound, upper_bound):
#     cut_array = array_to_cut[(lover_bound < boundary_array) & (boundary_array < upper_bound)]
#     return cut_array
#
#
# def inter_arm_cut(wl_arr = [] ,flux_arr = [], fluxerr_arr=[], transmission_arr=[], i_arr= [], start = [], end = []):
#
#     norm_length = 250
#     # Reformat
#     if i_arr == 'UVB':
#         low = 3100
#         high = 5550
#         wl_cut = cut(wl_arr, wl_arr, low, high)
#         flux_cut = cut(flux_arr, wl_arr, low, high)
#         fluxerr_cut = cut(fluxerr_arr, wl_arr, low, high)
#         transmission_cut = cut(transmission_arr, wl_arr, low, high)
#         end = np.median(flux_cut[-norm_length:])
#
#
#     if i_arr == 'VIS':
#         low = 5550
#         high = 10100
#         wl_cut = cut(wl_arr, wl_arr, low, high)
#         flux_cut = cut(flux_arr, wl_arr, low, high)
#         fluxerr_cut = cut(fluxerr_arr, wl_arr, low, high)
#         transmission_cut = cut(transmission_arr, wl_arr, low, high)
#         start = np.median(flux_cut[:norm_length])
#         flux_cut *= end / start
#         fluxerr_cut *= end / start
#         end = np.median(flux_cut[-norm_length:])
#
#     if i_arr == 'NIR':
#         low = 10100
#         high = 30700
#         wl_cut = cut(wl_arr, wl_arr, low, high)
#         flux_cut = cut(flux_arr, wl_arr, low, high)
#         fluxerr_cut = cut(fluxerr_arr, wl_arr, low, high)
#         transmission_cut = cut(transmission_arr, wl_arr, low, high)
#         start = np.median(flux_cut[:norm_length])
#         flux_cut *= end / start
#         fluxerr_cut *= end / start
#
#
#
#     return wl_cut, flux_cut, fluxerr_cut, transmission_cut, start, end



if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    import plotting.plotting as plot
    from xshoo.combine import inter_arm_cut

    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/')
    arms = ['UVB', 'VIS', 'NIR']
    redshifts = [1.1257, 1.98041, 1.57697, 1.82389, 1.51237, 2.096, 1.30914]
    for x,i in enumerate(sdssobjects):
        transmissionfiles = glob.glob(i+'*transmission*')
        obs = glob.glob(i+'*OBJECT*/*IDP*')
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
        for n in arms:
            print('In arm: '+n)
            obser = [k for k in obs if n in k]
            tran = [l for l in transmissionfiles if n in l]
            ob = fits.open(obser[0])

            tran = fits.open(tran[0])
            wl = 10.0*ob[1].data.field('WAVE')[0]



            flux = ob[1].data.field('FLUX')[0]
            err = ob[1].data.field('ERR')[0]
            transmission = tran[0].data

            wl, flux, err, transmission, start, end = inter_arm_cut(wl, flux, err, transmission, n, start, end)

            test.append(flux)


            fluxi = flux / transmission
            err = err / transmission

            wl_out.append(wl)
            flux_out.append(fluxi)
            err_out.append(err)

            if n == 'VIS':
                outfileop = obser[0]



        wl_out = np.hstack(wl_out)
        flux_out = np.hstack(flux_out)
        err_out = np.hstack(err_out)
        test = np.hstack(test)

        #Get SDSS spectrum
        imageSDSS= glob.glob(i+'*spec*.fits')
        data1SDSS= fits.open(imageSDSS[0])
        wl_sdss = 10.0**np.array(data1SDSS[1].data.field('loglam').flat)
        flux_sdss =  np.array(data1SDSS[1].data.field('flux').flat) * 1.e-17




        do_plot = True
        #Plotting
        if do_plot:
            fig, ax = pl.subplots()
            ax.plot(wl_out,flux_out, lw=0.1, color = 'black', label='Corrected')
            ax.plot(wl_out,test, lw=0.1, color = "red", label ="Uncorrected")
            ax.plot(wl_sdss,flux_sdss, lw=0.1, color = "blue", label ="SDSS")
            #Format axes
            #ax.set_xlim([3000,25000])
            #rng
            ax.set_xlim([4000,9000])
            #ax.set_ylim([-1e-15,4e-15])
            mask =  (wl_out > 4000) & (wl_out < 9000)

            ax.set_ylim([min(flux_out[mask] * 1.2),max(flux_out[mask] * 1.2)])
            ax.set_xlabel(r'Observed Wavelength  [\AA]')
            #ax.set_xlabel(r'Rest Wavelength  [\AA]')
            ax.set_ylabel(r'FLux [erg/cm$^2$/s/\AA]')
            ax.legend()
            composite = np.genfromtxt('linelist.txt', dtype=None)
            for p in range(len(composite)):
                xcoord = composite[p,][1]*(1+redshifts[x])
                mask = (wl_out > xcoord - 1) & (wl_out < xcoord + 1)
                y_val = np.median(flux_out[mask])
                ax.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.5)
                ax.annotate(composite[p,][0],xy=(xcoord, y_val * 1.2 ),fontsize='x-small')
            pl.title(obj_name)
            pl.tight_layout()
            file_name = "object"
            #pl.savefig(i+file_name+".pdf", clobber=True)
            pl.show(block=True)
#

        #Make equal length for saving
        wl_sdss = np.concatenate([wl_sdss,np.zeros(len(wl_out) - len(wl_sdss))])
        flux_sdss = np.concatenate([flux_sdss,np.zeros(len(flux_out) - len(flux_sdss))])

        #Saving to .dat file
        dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64), ("wl sdss", np.float64), ("flux sdss", np.float64) ]
        data = np.array(zip(wl_out, flux_out, err_out, wl_sdss, flux_sdss), dtype=dt)
        file_name = "Telluric_corrected_science"
        np.savetxt(i+file_name+".dat", data, header="wl flux fluxerror, sdss_wl, sdss_flux")#, fmt = ['%5.1f', '%2.15E'] )
     
            
            
            
            
            
            