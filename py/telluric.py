# -*- coding: utf-8 -*-





#np.set_printoptions(threshold=np.nan)



import ppxf
import ppxf_util as util

def find_best_template(wl, flux, err, hdr, spectral_library):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import splrep,splev
    from time import clock
    t2 = clock()
    # Read a galaxy spectrum and define the wavelength range
#    hdu = fits.open(telluric_observation)
    obs_spectrum = flux
    obs_spectrum_header = hdr

    obs_error_spectrum = err + 1
#    obs_badpix = hdu[2].data
#    print obs_badpix
    
#    wl_obs = 10.*(np.arange((np.shape(obs_spectrum)[0]))*obs_spectrum_header['CDELT1']+obs_spectrum_header['CRVAL1'])
    wl_obs = wl
#    obs_lambda_range = 10*(obs_spectrum_header['AWAVELMIN'] + np.array([0.,obs_spectrum_header['SPEC_BIN']*(obs_spectrum_header['NAXIS1']-1)]))
    obs_lambda_range = np.array([min(wl), max(wl)])
    # If the galaxy is at a significant redshift (z > 0.03), one would need to apply 
    # a large velocity shift in PPXF to match the template to the galaxy spectrum.
    # This would require a large initial value for the velocity (V > 1e4 km/s) 
    # in the input parameter START = [V,sig]. This can cause PPXF to stop! 
    # The solution consists of bringing the galaxy spectrum roughly to the 
    # rest-frame wavelength, before calling PPXF. In practice there is no 
    # need to modify the spectrum before the usual LOG_REBIN, given that a 
    # red shift corresponds to a linear shift of the log-rebinned spectrum. 
    # One just needs to compute the wavelength range in the rest-frame
    # and adjust the instrumental resolution of the galaxy observations.
    # This is done with the following three commented lines:
    
    z = 0.0 # Initial estimate of the galaxy redshift
    obs_lambda_range = obs_lambda_range/(1+z) # Compute approximate restframe wavelength range
   
    #logarithmically rebin while conserving flux    
    tell_obs, obs_lambda, velscale = util.log_rebin(obs_lambda_range, obs_spectrum)
    tell_obs_err, obs_lambda, velscale = util.log_rebin(obs_lambda_range, obs_error_spectrum)   
    
    #Normalize to avoid numerical issues
    tell_obs = tell_obs/np.median(tell_obs)
    tell_obs_err = tell_obs_err/np.median(tell_obs_err)

    #Load and prepare Model stellar library
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the target spectrum, to determine
    # the size needed for the array which will contain the template spectra.

    hdu = fits.open(spectral_library[0])
    library_spectrum = hdu[0].data
    library_spectrum_header = hdu[0].header
    wl_lib = np.e**(np.arange((library_spectrum_header['NAXIS1']))*library_spectrum_header['CDELT1']+library_spectrum_header['CRVAL1'])

    #Make empty template holder    
#    f = interpolate.interp1d(wl_lib,library_spectrum,kind='linear',bounds_error = False, fill_value=0.)
#    hdu_wave_short = f(wl_obs)
    
    f = splrep(wl_lib,library_spectrum, k=1)
    hdu_wave_short = splev(wl_obs,f)    
    
    lib_lambda_range = np.array([min(wl_obs),max(wl_obs)])
    tell_lib, lib_lambda, velscale = util.log_rebin(lib_lambda_range, hdu_wave_short, velscale=velscale)
    templates = np.empty((tell_lib.size,len(spectral_library)))

    # Convolve the whole library of spectral templates 
    # with the quadratic difference between the target and the 
    # library instrumental resolution. Logarithmically rebin 
    # and store each template as a column in the array TEMPLATES.
    
    for j in range(len(spectral_library)):
        t = clock() 
        hdu = fits.open(spectral_library[j])
        library_spectrum = hdu[0].data
        
        f = splrep(wl_lib,library_spectrum,k=1)
        interpolated_library_spectrum = splev(wl_obs,f)
        
        tell_lib, lib_lambda, velscale = util.log_rebin(lib_lambda_range,interpolated_library_spectrum, velscale=velscale)
        templates[:,j] = tell_lib/np.median(tell_lib) # Normalizes templates
        if j % 60 == 0:
            print 'Approximated remaining time (s) for setup of template spectra: '+ str(len(spectral_library)*(clock() - t) -  j*(clock() - t)) + 's' 

    #Excluding areas of strong telluric absorbtion from fitting
    if obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "UVB":        
        mask = (obs_lambda < np.log(5500)) & (obs_lambda > np.log(3100))
        goodPixels = np.where(mask == True)[0]
    elif obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "VIS":
        mask = (obs_lambda > np.log(5500)) & (obs_lambda < np.log(6860)) | (obs_lambda > np.log(7045)) & (obs_lambda < np.log(7140)) | (obs_lambda > np.log(7355)) & (obs_lambda < np.log(7570)) | (obs_lambda > np.log(7710)) & (obs_lambda < np.log(8090)) | (obs_lambda > np.log(8400)) & (obs_lambda < np.log(8900)) | (obs_lambda > np.log(9900)) & (obs_lambda < np.log(10100))
        goodPixels = np.where(mask == True)[0]     
    elif obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "NIR":   
        mask = (obs_lambda < np.log(10950)) | (obs_lambda > np.log(12240)) & (obs_lambda < np.log(12500)) | (obs_lambda > np.log(12800)) & (obs_lambda < np.log(12950)) | (obs_lambda > np.log(15300)) & (obs_lambda < np.log(17100)) | (obs_lambda > np.log(21000)) & (obs_lambda < np.log(23700))  
        goodPixels = np.where(mask == True)[0]
        
        
    # The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV. This assume the redshift is negligible.
    # In the case of a high-redshift galaxy one should de-redshift its 
    # wavelength to the rest frame before using the line below (see above).

    c = 299792.458
    dv = 0# (logLam2[0]-logLam1[0])*c # km/s
    vel = 0 # Initial estimate of the galaxy velocity in km/s
    start = [vel, 50.] # (km/s), starting guess for [V,sigma]
    
    # Here the actual fit starts. The best fit is plotted on the screen.
    # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
    pp = ppxf.ppxf(templates, tell_obs, tell_obs_err, velscale, start,
                       goodpixels=goodPixels, plot=False, moments=4, 
                       degree=4, mdegree=4, vsyst=dv, clean=True, regul=1)
#    plt.show(block=True)
    print "Formal errors:"    
    print "     dV    dsigma   dh3      dh4"
    print "".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2))
    
    print 'elapsed time (s) for ppxf: ', clock() - t2   
    
    # If the galaxy is at significant redshift z and the wavelength has been
    # de-redshifted with the three lines "z = 1.23..." near the beginning of 
    # this procedure, the best-fitting redshift is now given by the following 
    # commented line (equation 2 of Cappellari et al. 2009, ApJ, 704, L34):
    
    print 'Best-fitting redshift z:', (z + 1)*(1 + pp.sol/c) - 1
    print 'Best-fitting error on redshift z:',((z + 1)*(1 + pp.sol/c) - 1) - ((z + 1)*(1 + pp.error*np.sqrt(pp.chi2)/c) - 1) 
    print 'Rebinning to linear axes'
    import spec   
    obj_spec = spec.resamplespec(wl_obs,np.e**obs_lambda,pp.galaxy, oversamp =1000)
    template_fit = spec.resamplespec(wl_obs,np.e**obs_lambda,pp.bestfit, oversamp =1000)
    plt.plot(wl_obs,obj_spec, color='black', linestyle = 'steps-mid')
    plt.plot(wl_obs, template_fit, color='red')
    return obj_spec,template_fit,obs_spectrum_header
#------------------------------------------------------------------------------



# def inter_arm_cut(wl_arr = [] ,flux_arr = [], fluxerr_arr=[], i_arr= []):
#     wl = 0
#     flux = 0
#     fluxerr = 0
#
#     #Reformat
#     if i_arr == 'UVB':
#         upper = 5550
#         lower = 3100
#         wl = wl_arr[(lower < wl_arr) & (wl_arr < upper)]
#         flux = flux_arr[(lower < wl_arr) & (wl_arr < upper)]
#         fluxerr = fluxerr_arr[(lower < wl_arr) & (wl_arr < upper)]
#
#     if i_arr == 'VIS':
#         upper = 10100
#         lower = 5550
#         wl = wl_arr[(lower < wl_arr) & (wl_arr < upper)]
#         flux = flux_arr[(lower < wl_arr) & (wl_arr < upper)]
#         fluxerr = fluxerr_arr[(lower < wl_arr) & (wl_arr < upper)]
#
#     if i_arr == 'NIR':
#         lower = 10100
#         wl = wl_arr[(lower< wl_arr)]
#         flux = flux_arr[(lower < wl_arr)]
#         fluxerr = fluxerr_arr[(lower < wl_arr)]
#
#     return wl, flux, fluxerr



    
if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobject = 'SDSS1437-0147'
    sdssobjects = glob.glob(root_dir+'*/')
    #Load in Model steller spetra
    #library = glob.glob('/Users/jselsing/nosync/spectral_libraries/phoenix_spectral_library/TRIAL/*.fits')
    arms = ['UVB', 'VIS', 'NIR']
    #arms = ['NIR']
    library = glob.glob('/Users/jselsing/nosync/spectral_libraries/phoenix_spectral_library/R10000RES/*/*.fits')


    for i in sdssobjects:
        print 'Working on object: '+i
        tell_files = glob.glob(i+'TELLURIC_STAR/*IDP*')
        
        
        for n in arms:
            print 'In arm: '+n
            tell_file = [k for k in tell_files if n in k]
            tell_file = fits.open(tell_file[0])
            wl = 10.0*tell_file[1].data.field('WAVE')[0]
            flux = tell_file[1].data.field('FLUX')[0]
            err = tell_file[1].data.field('ERR')[0]
            # wl, flux, err = inter_arm_cut(wl_arr=wl, flux_arr=flux,
            #                                 fluxerr_arr=err, i_arr=n)


            gal, fit, hdr = find_best_template(wl, flux, err, tell_file[0].header, library)
            trans = (gal/fit)
            fileout = np.array([wl, trans])
            print(fileout)

            fits.writeto(i+'transmission_'+n+'.fits',trans, hdr, clobber=True)

            dt = [("wl", np.float64), ("telluric_star", np.float64), ("Optimal_template_fit", np.float64)]
            data = np.array(zip(wl, gal, fit), dtype=dt)
            file_name = "Telluric_correction_QC"
            np.savetxt(i+file_name+n+".dat", data, header="wl telluric_star Optimal_template_fit")#, fmt = ['%5.1f', '%2.15E'] )            
            

            # print "close the plot to continue"
            # pl.show(block=True)
    

    