# -*- coding: utf-8 -*-


import ppxf
import ppxf_util as util

def find_best_template(wl_obs, flux, err, hdr, spectral_library):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import splrep,splev
    from time import clock
    t2 = clock()

    # Read a spectrum and define the wavelength range
    obs_spectrum = flux
    obs_spectrum_header = hdr
    obs_error_spectrum = err

    obs_lambda_range = np.array([min(wl), max(wl)])
    z = 0.0 # Initial estimate of the galaxy redshift
    obs_lambda_range = obs_lambda_range/(1+z) # Compute approximate restframe wavelength range

    #Get median of positive values
    m = np.median(obs_error_spectrum[obs_error_spectrum > 0])
    # Assign the median to the negative elements
    obs_error_spectrum[obs_error_spectrum <= 0] = m




    #logarithmically rebin while conserving flux    
    tell_obs, obs_lambda, velscale = util.log_rebin(obs_lambda_range, obs_spectrum)
    tell_obs_err, obs_lambda, velscale = util.log_rebin(obs_lambda_range, obs_error_spectrum)   
    
    #Normalize to avoid numerical issues
    norm = np.median(tell_obs)
    tell_obs = tell_obs/norm
    tell_obs_err = tell_obs_err/norm

    # Load and prepare Model stellar library
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the target spectrum, to determine
    # the size needed for the array which will contain the template spectra.

    hdu = fits.open(spectral_library[0])
    library_spectrum = hdu[0].data
    library_spectrum_header = hdu[0].header
    wl_lib = np.e**(np.arange((library_spectrum_header['NAXIS1']))*library_spectrum_header['CDELT1']+library_spectrum_header['CRVAL1'])

    #Make empty template holder
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

        # Interpolate template spectrum to match input spectrum
        f = splrep(wl_lib,library_spectrum,k=1)
        interpolated_library_spectrum = splev(wl_obs,f)

        # Logarithmically rebin template spectra
        tell_lib, lib_lambda, velscale = util.log_rebin(lib_lambda_range,interpolated_library_spectrum, velscale=velscale)
        templates[:,j] = tell_lib/np.median(tell_lib) # Normalizes templates
        if j % 60 == 0:
            print 'Approximated remaining time (s) for setup of template spectra: '+ str(len(spectral_library)*(clock() - t) -  j*(clock() - t)) + 's' 

    #Excluding areas of strong telluric absorbtion from fitting
    if obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "UVB":        
        mask = (obs_lambda < np.log(5500)) & (obs_lambda > np.log(3100))
        goodPixels = np.where(mask == True)[0]
    elif obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "VIS":
        mask = (obs_lambda > np.log(5500)) & (obs_lambda < np.log(6350)) | (obs_lambda > np.log(6380)) & (obs_lambda < np.log(6860)) | (obs_lambda > np.log(7045)) & (obs_lambda < np.log(7140)) | (obs_lambda > np.log(7355)) & (obs_lambda < np.log(7570)) | (obs_lambda > np.log(7710)) & (obs_lambda < np.log(8090)) | (obs_lambda > np.log(8400)) & (obs_lambda < np.log(8900)) | (obs_lambda > np.log(9900)) & (obs_lambda < np.log(10100))
        goodPixels = np.where(mask == True)[0]     
    elif obs_spectrum_header['HIERARCH ESO SEQ ARM'] == "NIR":   
        mask = (obs_lambda > np.log(10000)) & (obs_lambda < np.log(10950)) | (obs_lambda > np.log(12240)) & (obs_lambda < np.log(12500)) | (obs_lambda > np.log(12800)) & (obs_lambda < np.log(12950)) | (obs_lambda > np.log(15300)) & (obs_lambda < np.log(17100)) | (obs_lambda > np.log(21000)) & (obs_lambda < np.log(23700))
        goodPixels = np.where(mask == True)[0]
        

    # Initial parameters for LOSVD
    c = 299792.458
    dv = 0# (logLam2[0]-logLam1[0])*c # km/s
    vel = 0 # Initial estimate of the LOSVD km/s
    start = [vel, 10.] # (km/s), starting guess for [V,sigma]
    
    # Here the actual fit starts.
    pp = ppxf.ppxf(templates, tell_obs, tell_obs_err, velscale, start,
                       goodpixels=goodPixels, plot=True, moments=4,
                       degree=5, mdegree=5, vsyst=dv, regul = 10)

    print "Formal errors:"    
    print "     dV    dsigma   dh3      dh4"
    print "".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2))
    
    print 'elapsed time (s) for ppxf: ', clock() - t2   
    

    
    print 'Best-fitting redshift z:', (z + 1)*(1 + pp.sol/c) - 1
    print 'Best-fitting error on redshift z:',((z + 1)*(1 + pp.sol/c) - 1) - ((z + 1)*(1 + pp.error*np.sqrt(pp.chi2)/c) - 1)


    print 'Rebinning to linear axes'
    import spec   
    obj_spec = spec.resamplespec(wl_obs,np.e**obs_lambda,pp.galaxy, oversamp = 1000)
    template_fit = spec.resamplespec(wl_obs,np.e**obs_lambda,pp.bestfit, oversamp =1000)
    # plt.plot(wl_obs,obj_spec, color='black', linestyle = 'steps-mid')
    # plt.plot(wl_obs, template_fit, color='red')

    return obj_spec,template_fit,obs_spectrum_header
#------------------------------------------------------------------------------



    
if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np
    from scipy.interpolate import splrep,splev
    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobject = glob.glob(root_dir+'SDSS1431+0535/')
    sdssobjects = glob.glob(root_dir+'*/')
    # sdssobjects = ['SDSS0209-0947/', 'SDSS0043+0114/', 'SDSS0155-1023/',]
    # sdssobjects = [ 'SDSS0303+0027/', 'SDSS0323-0029/', 'SDSS0842+0151/', 'SDSS1002+0331/' , 'SDSS1158-0322/']
    # sdssobjects = [ 'SDSS0820+1306/', 'SDSS1150-0023/', 'SDSS1219-0100/', 'SDSS1236-0331/' , 'SDSS1354-0013/',
    #                 'SDSS1431+0535/', 'SDSS1437-0147/']
    sdssobjects = ['SDSS1150-0023/']
    #Load in Model steller spetra

    # arms = [ 'VIS', 'NIR']
    arms = ['NIR']

    library = glob.glob('/Users/jselsing/nosync/spectral_libraries/phoenix_spectral_library/R10000RES/*/*.fits')
    # library = glob.glob('/Users/jselsing/nosync/spectral_libraries/phoenix_spectral_library/TEMPLATES/*/*.fits')

    for i in sdssobjects:
        i = root_dir + i
        print 'Working on object: '+i
        tell_files = glob.glob(i+'TELLURIC_STAR/*IDP*')
        # print(tell_files)
        # master_response = glob.glob(i+'M.X*.fits')
        
        for l,n in enumerate(arms):
            print 'In arm: '+n
            tell_file = [k for k in tell_files if n in k]
            tell_file = fits.open(tell_file[0])
            wl = 10.0*tell_file[1].data.field('WAVE')[0]
            # resp = fits.open(master_response[l+1])
            # response_wl = resp[1].data.field('LAMBDA')*10.0
            # response = (resp[1].data.field('RESPONSE'))
            # interp = splrep(response_wl,response)
            # response_interp = splev(wl,interp)





            flux = tell_file[1].data.field('FLUX')[0] #* response_interp
            err = tell_file[1].data.field('ERR')[0] #* response_interp

            gal, fit, hdr = find_best_template(wl, flux, err, tell_file[0].header, library)
            trans = (gal/fit)


            fits.writeto(i+'transmission_'+n+'.fits',trans, hdr, clobber=True)

            dt = [("wl", np.float64), ("telluric_star", np.float64), ("Optimal_template_fit", np.float64)]
            data = np.array(zip(wl, gal, fit), dtype=dt)
            file_name = "Telluric_correction_QC"
            np.savetxt(i+file_name+n+".dat", data, header="wl telluric_star Optimal_template_fit")
            

            print("close the plot to continue")
            pl.show(block=True)
    

    