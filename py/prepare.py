# -*- coding: utf-8 -*-
"""
Call as Analyser.analysis(ebv,z,sdss,name)
"""
#from astropy.io import fits
#import glob
#import numpy as np
#import matplotlib.pyplot as plt
#from gen_methods import smooth
#
#
#from scipy import interpolate
#from scipy import ndimage





#def prepare(sdssobject)




def analysis(wl, flux, err, transmission, ebv, init_z, sdss, name, tell_corr_method='synthetic'):
    from unred import ccm_unred


    arm = ['UVB','VIS','NIR']
    for n in arm:
     
        #Deredden using E(B-V)-value from Schlegel map
        flux = ccm_unred(wl, flux , ebv)
        
        #Load transmission curve generated from synthetic steller spectra
#        trans = fits.open(name+'/sub'+n+'.fits')
        
        #Load transmission curve generated from direct method
#        transdir = fits.open(name+'/sub'+n+'dir.fits')
        

        
        
        #Masking bad values
#        mask = (fitsfile[0].data > 0) & (trans[0].data > 0) & (transdir[0].data > 0 ) & (fitsfile[0].data - fitsfile[1].data > 0)   

            
        if tell_corr_method == 'synth':
            if n == 'UVB':
                mask = wl < 5575
                wlUVB = wl[mask]
                wlUVB /= (1.0 + z)
                fluxUVB = (fitsfile[0].data)[mask]
                fluxerrUVB = (fitsfile[0].data)[mask]
                flux_corrUVB = (fitsfile[0].data/trans[0].data)[mask]
                flux_correrrUVB = (fitsfile[1].data/trans[0].data)[mask]
                
            if n == 'VIS':
                mask = (wl >= 5575) & (wl < 10100)
                wlVIS = wl[mask]
                wlVIS /= (1.0 + z)
                fluxVIS = (fitsfile[0].data)[mask]
                fluxerrVIS = (fitsfile[0].data)[mask]
                flux_corrVIS = (fitsfile[0].data/trans[0].data)[mask]
                flux_correrrVIS = (fitsfile[1].data/trans[0].data)[mask]
        
            if n == 'NIR':
                mask = wl >= 10100
                wlNIR = wl[mask]
                wlNIR /= (1.0 + z)
                fluxNIR = (fitsfile[0].data)[mask]
                fluxerrNIR = (fitsfile[0].data)[mask]
                flux_corrNIR = (fitsfile[0].data/trans[0].data)[mask]
                flux_correrrNIR = (fitsfile[1].data/trans[0].data)[mask]
        
        if tell_corr_method == 'direct':
            if n == 'UVB':
                mask = wl < 5575
                wlUVB = wl[mask]
                wlUVB /= (1.0 + z)
                fluxUVB = (fitsfile[0].data)[mask]
                fluxerrUVB = (fitsfile[0].data)[mask]
                flux_corrUVB = (fitsfile[0].data/transdir[0].data)[mask]
                flux_correrrUVB = (fitsfile[1].data/transdir[0].data)[mask]
                
            if n == 'VIS':
                mask = (wl >= 5575) & (wl < 10100)
                wlVIS = wl[mask]
                wlVIS /= (1.0 + z)
                fluxVIS = (fitsfile[0].data)[mask]
                fluxerrVIS = (fitsfile[0].data)[mask]
                flux_corrVIS = (fitsfile[0].data/transdir[0].data)[mask]
                flux_correrrVIS = (fitsfile[1].data/transdir[0].data)[mask]
        
            if n == 'NIR':
                mask = wl >= 10100
                wlNIR = wl[mask]
                wlNIR /= (1.0 + z)
                fluxNIR = (fitsfile[0].data)[mask]
                fluxerrNIR = (fitsfile[0].data)[mask]
                flux_corrNIR = (fitsfile[0].data/transdir[0].data)[mask]
                flux_correrrNIR = (fitsfile[1].data/transdir[0].data)[mask]

              
    from itertools import chain
    wl = np.asarray(list(chain(wlUVB,wlVIS,wlNIR)))
    flux = np.asarray(list(chain(fluxUVB,fluxVIS,fluxNIR)))
    fluxerr = np.asarray(list(chain(fluxerrUVB,fluxerrVIS,fluxerrNIR)))
    flux_corr = np.asarray(list(chain(flux_corrUVB,flux_corrVIS,flux_corrNIR)))
    flux_correrr = np.asarray(list(chain(flux_correrrUVB,flux_correrrVIS,flux_correrrNIR)))
    
    return wl, flux , fluxerr, flux_corr, flux_correrr

def inter_arm_cut(wl_arr = [] ,flux_arr = [], fluxerr_arr=[], transmission_arr=[], i_arr= []):
    wl = 0
    flux = 0
    fluxerr = 0
    transmission = 0

    #Reformat
    if i_arr == 'UVB':
        wl = wl_arr[(3100 < wl_arr) & (wl_arr < 5550)]
        flux = flux_arr[(3100 < wl_arr) & (wl_arr < 5550)]
        fluxerr = fluxerr_arr[(3100 < wl_arr) & (wl_arr < 5550)]
        transmission = transmission_arr[(3100 < wl_arr) & (wl_arr < 5550)]
    
    if i_arr == 'VIS':
        wl = wl_arr[(5550 < wl_arr) & (wl_arr < 10100)]
        flux = flux_arr[(5550 < wl_arr) & (wl_arr < 10100)]
        fluxerr = fluxerr_arr[(5550 < wl_arr) & (wl_arr < 10100)]
        transmission = transmission_arr[(5550 < wl_arr) & (wl_arr < 10100)]
        
    if i_arr == 'NIR':
        upper = 30700
        lower = 10100
        wl = wl_arr[(lower< wl_arr) & (wl_arr < upper)]
        flux = flux_arr[(lower < wl_arr) & (wl_arr < upper)]
        fluxerr = fluxerr_arr[(lower < wl_arr) & (wl_arr < upper)]
        transmission = transmission_arr[(lower < wl_arr) & (wl_arr < upper)]
        
    return wl, flux, fluxerr, transmission





    
if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    import plotting.plotting as plot
    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*/')
    arms = ['UVB', 'VIS', 'NIR']     
    for i in sdssobjects:
        transmissionfiles = glob.glob(i+'*tran*')
        obs = glob.glob(i+'*OBJECT*/*IDP*')
        obj_name = i[-14:-1]

        wl_out = []
        flux_out = []
        test = []
        err_out = []
        for n in arms:
#            print 'In arm: '+n
            ob = [k for k in obs if n in k]
            tran = [l for l in transmissionfiles if n in l]

            ob = fits.open(ob[0])

            tran = fits.open(tran[0])

            wl = 10.0*ob[1].data.field('WAVE')[0]
            flux = ob[1].data.field('FLUX')[0]
            err = ob[1].data.field('ERR')[0]
            transmission = tran[0].data
            wl, flux, err, transmission = inter_arm_cut(wl_arr=wl, flux_arr=flux,
                                                        fluxerr_arr=err, transmission_arr=transmission, i_arr=n)
                                                        
            test.append(flux)
            
            fluxi = flux / transmission
            err /= transmission

            wl_out.append(wl)
            flux_out.append(fluxi)
            err_out.append(err)
            

            
        wl_out = np.hstack(wl_out)
        flux_out = np.hstack(flux_out)       
        err_out = np.hstack(err_out)   
        test = np.hstack(test)
        fig = plot.plot_data(wl_out,flux_out,xrng=[3000,25000], yrng=[-1e-15,4e-15], title = str(obj_name), lw=0.2)
        fig = plot.plot_data(wl_out,test,xrng=[3000,25000], yrng=[-1e-15,4e-15], title = str(obj_name), lw=0.2, color = "red", fig=fig)
#        pl.plot(wl_out,flux_out, lw=0.5, color="black")
#        pl.plot(wl_out,test, lw=0.5, color="red")
        
        pl.show(block=True)
  
        dt = [("wl", np.float64), ("flux", np.float64), ("error", np.float64)]
        data = np.array(zip(wl_out, flux_out, err_out), dtype=dt)
        file_name = "Telluric_corrected_science"
        np.savetxt(i+file_name+".dat", data, header="wl flux fluxerror")#, fmt = ['%5.1f', '%2.15E'] )
     
            
            
            
            
            
            