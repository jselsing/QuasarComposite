# -*- coding: utf-8 -*-
"""
Call as Analyser.analysis(ebv,z,sdss,name)
"""
from __future__ import division, print_function

__all__ = []
__author__ = 'jselsing'

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')


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
    sdssobjects = ['SDSS0021+0043/', 'SDSS0022+0124/', 'SDSS0820+1306/', 'SDSS1013+0851/', 'SDSS1101+0548/', 'SDSS1150-0023/', 'SDSS1219-0100/', 'SDSS1236-0331/', 'SDSS1249-0559/', 'SDSS1354-0013/', 'SDSS1431+0535/', 'SDSS1437-0147/', 'SDSS2123-0050/', 'SDSS2313+0034/']
    sdssobjects = ['SDSS0820+1306/', 'SDSS1150-0023/', 'SDSS1219-0100/', 'SDSS1236-0331/', 'SDSS1354-0013/', 'SDSS1431+0535/', 'SDSS1437-0147/']
    sdssobjects = ['SDSS1431+0535/']

    arms = ['UVB', 'VIS', 'NIR']
    redshifts = [1.1257, 1.98041, 1.57697, 1.82389, 1.51237, 2.096, 1.30914]
    for x,i in enumerate(sdssobjects):
        transmissionfiles = glob.glob(root_dir+i+'*transmission*')
        obs = glob.glob(root_dir+i+'*OBJECT*/*IDP*')
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
            # pl.plot(wl,flux)

            pl.plot(wl,flux, lw = 0.5, linestyle='steps-mid')


            #trans_new = []
            for j, k in enumerate(transmission):
                if k <= 1e-10:
                    transmission[j] = 1
            print(min(transmission))
            #fluxerr = np.array(fluxerr_new)

            fluxi = flux / transmission
            err = err / transmission
            pl.plot(wl,transmission*np.median(flux), lw = 0.5, linestyle='steps-mid')
            pl.plot(wl,fluxi, lw = 0.5, linestyle='steps-mid')
            pl.show()

            wl_tmp, fluxi, err_tmp, transmission_tmp, start_tmp, end_tmp = inter_arm_cut(wl, fluxi, err, transmission, n, start, end)
            wl, flux, err, transmission, start, end = inter_arm_cut(wl, flux, err, transmission, n, start, end)
            test.append(flux)




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
        imageSDSS= glob.glob(root_dir+i+'*spec*.fits')
        data1SDSS= fits.open(imageSDSS[0])
        wl_sdss = 10.0**np.array(data1SDSS[1].data.field('loglam').flat)
        flux_sdss =  np.array(data1SDSS[1].data.field('flux').flat) * 1.e-17




        do_plot = False
        #Plotting
        if do_plot:
            fig, ax = pl.subplots()
            ax.plot(wl_out,flux_out, lw=0.1, color = 'black', label='Corrected')
            ax.plot(wl_out,test, lw=0.1, color = "red", label ="Uncorrected")
            # ax.plot(wl_sdss,flux_sdss, lw=0.1, color = "blue", label ="SDSS")
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
            # composite = np.genfromtxt('linelist.txt', dtype=None)
            # for p in range(len(composite)):
            #     xcoord = composite[p,][1]*(1+redshifts[x])
            #     mask = (wl_out > xcoord - 1) & (wl_out < xcoord + 1)
            #     y_val = np.median(flux_out[mask])
            #     ax.axvline(x=xcoord,color='green',linestyle='dashed', lw=0.5)
            #     ax.annotate(composite[p,][0],xy=(xcoord, y_val * 1.2 ),fontsize='x-small')
            pl.title(obj_name)
            pl.tight_layout()
            file_name = "object"
            #pl.savefig(i+file_name+".pdf", clobber=True)
            pl.show(block=True)
#

        #Make equal length for saving
        wl_sdss = np.concatenate([wl_sdss,np.zeros(len(wl_out) - len(wl_sdss))])
        flux_sdss = np.concatenate([flux_sdss,np.zeros(len(flux_out) - len(flux_sdss))])

        # fluxerr_new = []
        # for j, (k, l) in enumerate(zip(flux_out,err_out)):
        #     if k > 2 * flux_out[j-2] and k > 0:
        #         fluxerr_new.append(100*l)
        #     elif k < 1/2 * flux_out[j-2] and k > 0:
        #         fluxerr_new.append(100*l)
        #     else:
        #         fluxerr_new.append(l)
        # from gen_methods import smooth
        # err_out = smooth(np.array(fluxerr_new), window_len=5, window='hanning')


        # fluxerr_new = []
        # for j , (k, l) in enumerate(zip(flux_out[:-1],err_out[:-1])):
        #     if k > 1.5 * flux_out[j-1] and k > 0:
        #         fluxerr_new.append(abs(1000*err_out[j+1]))
        #     elif k < 0.75 * flux_out[j-1] and k > 0:
        #         fluxerr_new.append(abs(1000*err_out[j+1]))
        #     else:
        #         fluxerr_new.append(abs(err_out[j+1]))
        # fluxerr_new.append(0)
        # from gen_methods import smooth
        # err_out = smooth(np.array(fluxerr_new), window_len=5, window='hanning')


        # fluxerr_new = []
        # for j , (k, l) in enumerate(zip(flux_out[:-1],err_out[:-1])):
        #     if k > 1.2 * flux_out[j-1] or k < 0:
        #        fluxerr_new.append(abs(1000*err_out[j+1]))
        #     elif k < 0.80 * flux_out[j-1] or k < 0:
        #         fluxerr_new.append(abs(1000*err_out[j+1]))
        #     else:
        #         fluxerr_new.append(abs(err_out[j+1]))
        # fluxerr_new.append(0)
        # from gen_methods import smooth
        # fluxerr_new = smooth(np.array(fluxerr_new), window_len=1, window='hanning')


        bp_map = []
        for j , (k, l) in enumerate(zip(flux_out[:-1],err_out[:-1])):
            if k > 1.1 * flux_out[j-1] or k < 0:
               bp_map.append(1)
            elif k < 0.90 * flux_out[j-1] or k < 0:
               bp_map.append(1)
            else:
               bp_map.append(0)
        bp_map.append(1)






        # print(len(wl_out))
        # print(len(fluxerr_new))
        # print(len(flux_out))
        # pl.plot(wl_out, flux_out, linestyle='steps-mid')
        # pl.plot(wl_out, bp_map, linestyle='steps-mid')
        # pl.show()


        #Saving to .dat file
        dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64), ("bp map", np.float64), ("wl sdss", np.float64), ("flux sdss", np.float64) ]
        data = np.array(zip(wl_out, flux_out, err_out, bp_map, wl_sdss, flux_sdss), dtype=dt)
        file_name = "Telluric_corrected_science"
        np.savetxt(root_dir+i+file_name+".dat", data, header="wl flux fluxerror bp_map sdss_wl, sdss_flux")#, fmt = ['%5.1f', '%2.15E'] )
     
            
            
            
            
            
            