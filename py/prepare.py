# -*- coding: utf-8 -*-
"""
Call as Analyser.analysis(ebv,z,sdss,name)
"""
from astropy.io import fits
import glob
import numpy as np
import matplotlib.pyplot as plt
from gen_methods import smooth


from scipy import interpolate
from scipy import ndimage

from unred import ccm_unred




def analysis(ebv,z,sdss,name, printsn = False, plot=False,color = None, errorbars = False,offset = 1,tell_corr_method='direct',ax=None):
    arm = ['UVB','VIS','NIR']
    for n in arm:
        #Reading in the reduction products
        file_path =glob.glob(name+'/OBJECT/*FLUX*'+n+'*.fits')[0]
        fitsfile = fits.open(file_path)
        
        #Reading in telluric star
        synth_path = glob.glob(name+'/TELLURIC_STAR/*1D*'+n+'*.fits')[0]
        fitssynth = fits.open(synth_path)
        resp = fits.open(name+'/RESPONSE_MERGE1D_SLIT_'+n+'.fits')[1].data.field('RESPONSE')
        fitssynth[0].data *= resp
        
        #Make wavelength array
        wl = 10.*(np.arange((np.shape(fitsfile[0].data)[0]))*fitsfile[0].header['CDELT1']+fitsfile[0].header['CRVAL1'])
     
        #Deredden using E(B-V)-value from Schlegel map
        fitsfile[0].data = ccm_unred(wl, fitsfile[0].data , ebv)
        
        #Load transmission curve generated from synthetic steller spectra
        trans = fits.open(name+'/sub'+n+'.fits')
        
        #Load transmission curve generated from direct method
        transdir = fits.open(name+'/sub'+n+'dir.fits')
        
        #synthetic transmission curves
    #    skymodel = fits.open('skytableR7450.fits')
    #    lammodel = skymodel[1].data.field('lam') * 1e4
    #    transmodel = skymodel[1].data.field('trans')
    #    transmodel = ndimage.gaussian_filter1d(transmodel, sigma= 0.5)
    #    f = interpolate.UnivariateSpline(lammodel,transmodel,k=1,s=0)
    #    if n== 'UVB':
    #        transmodel = f(wl*(1+77.3/3e5))
    #    elif n == 'VIS':     
    #        transmodel = f(wl*(1+76.3/3e5))
    #    elif n == 'NIR':    
    #        transmodel = f(wl*(1+96.3/3e5))
    #    else:
    #        print 'Strange things are amiss'
        
        
        
        
        
        
        bins= 5
#        fitsfile[0].data = binning1d(fitsfile[0].data,bins) 
#        fitsfile[1].data = binning1d(fitsfile[1].data,bins)  
#        trans[0].data = binning1d(trans[0].data,bins)  
#        transdir[0].data = binning1d(transdir[0].data,bins)     
#        wl = binning1d(wl,bins)         
        
        
        
        
        
        
        #Masking bad values
        mask = (fitsfile[0].data > 0) & (trans[0].data > 0) & (transdir[0].data > 0 ) & (fitsfile[0].data - fitsfile[1].data > 0)   
            #Excluding areas of strong telluric absorbtion from fitting
#        if n == "UVB":        
#            mask = (wl < np.log(5575)) & (wl > np.log(3500)) 
#    
#        if n == "VIS":
#            mask = (wl > np.log(5500)) & (wl < np.log(6860)) | (wl > np.log(7045)) & (wl < np.log(7140)) | (wl > np.log(7355)) & (wl < np.log(7570)) | (wl > np.log(7710)) & (wl < np.log(8090)) | (wl > np.log(8400)) & (wl < np.log(8900)) | (wl > np.log(9900)) & (wl < np.log(10100))
#     
#        if n == "NIR":   
#            mask = (wl < np.log(10950)) | (wl > np.log(12240)) & (wl < np.log(12500)) | (wl > np.log(12800)) & (wl < np.log(12950)) | (wl > np.log(15300)) & (wl < np.log(17100)) | (wl > np.log(21000)) & (wl < np.log(21700))  

        
        wl =  wl[mask]
        trans[0].data = trans[0].data[mask]
        transdir[0].data = transdir[0].data[mask]
    #    transmodel = transmodel[mask]
        fitsfile[0].data = fitsfile[0].data[mask]
        fitsfile[1].data = fitsfile[1].data[mask]
       
        #Calculate and prin S/N for various wavelength regions
        if printsn == True:
            if n == 'UVB':
                print 'S/N for various regions of heavy telluric absorption in '+name
            if n== 'VIS':
                print '(wlVIS > 7500) & (wlVIS < 7750)'
                mask = (wl > 7500) & (wl < 7750)
                
                print np.mean((fitsfile[0].data)[mask])/np.std((fitsfile[0].data)[mask])
                print np.mean((fitsfile[0].data/trans[0].data)[mask])/np.std((fitsfile[0].data/trans[0].data)[mask])
                print np.mean((fitsfile[0].data/transdir[0].data)[mask])/np.std((fitsfile[0].data/transdir[0].data)[mask])
                
                print '(wl > 9250) & (wl < 9750)'
                mask = (wl > 9250) & (wl < 9750)
                
                print np.mean((fitsfile[0].data)[mask])/np.std((fitsfile[0].data)[mask])
                print np.mean((fitsfile[0].data/trans[0].data)[mask])/np.std((fitsfile[0].data/trans[0].data)[mask])
                print np.mean((fitsfile[0].data/transdir[0].data)[mask])/np.std((fitsfile[0].data/transdir[0].data)[mask])
                
            if n == 'NIR':        
                
                
                print '(wl > 10750) & (wl < 11750)'
                mask = (wl > 10750) & (wl < 11750)
                
                print np.mean((fitsfile[0].data)[mask])/np.std((fitsfile[0].data)[mask])
                print np.mean((fitsfile[0].data/trans[0].data)[mask])/np.std((fitsfile[0].data/trans[0].data)[mask])
                print np.mean((fitsfile[0].data/transdir[0].data)[mask])/np.std((fitsfile[0].data/transdir[0].data)[mask])
                
                print '(wl > 12500) & (wl < 15500)'
                mask = (wl > 12500) & (wl < 15500)
                
                print np.mean((fitsfile[0].data)[mask])/np.std((fitsfile[0].data)[mask])
                print np.mean((fitsfile[0].data/trans[0].data)[mask])/np.std((fitsfile[0].data/trans[0].data)[mask])
                print np.mean((fitsfile[0].data/transdir[0].data)[mask])/np.std((fitsfile[0].data/transdir[0].data)[mask])
                
                print '(wl > 19950) & (wl < 20250)'
                mask = (wl > 19950) & (wl < 20250)
                
                print np.mean((fitsfile[0].data)[mask])/np.std((fitsfile[0].data)[mask])
                print np.mean((fitsfile[0].data/trans[0].data)[mask])/np.std((fitsfile[0].data/trans[0].data)[mask])
                print np.mean((fitsfile[0].data/transdir[0].data)[mask])/np.std((fitsfile[0].data/transdir[0].data)[mask])
            
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
            
        #Redshifting to object restframe    
        wl /= (1.0 + z)
        
       
        
        
        #plotting
        if plot == True:
#            fig, ax = plt.subplots()
#            ax.plot(wl,offset*fitsfile[0].data, drawstyle='steps-mid',color = color,lw=0.3)
            if tell_corr_method == 'synth':
                ax.plot(wl,offset*fitsfile[0].data/trans[0].data,drawstyle='steps-mid',color = color,lw=0.3)
            if tell_corr_method == 'direct':
                ax.plot(wl,offset*fitsfile[0].data/transdir[0].data,drawstyle='steps-mid',color = color,lw=0.3)

            #ax.step(wl,fitsfile[0].data/transmodel, color='seagreen',alpha=1.0,label='Corrected '+n+'model')


             #Overplot SDSS#
            imageSDSS= name+'/'+sdss
            data1SDSS,hdr1SDSS = fits.getdata(imageSDSS, header=True)
            wl = 10**np.array(data1SDSS.field('loglam').flat)/(1.+z)
            flux =  np.array(data1SDSS.field('flux').flat) * 1.e-17
            ax.step(wl,offset*smooth(flux, window_len=1, window='hanning'),'m-')#,label='SDSS')
            
             #Overplot Johans composite
#            vandenberk = np.loadtxt(name+'/composite.txt')
#            mask = (wl > 1600) & (wl < 1700) 
#            maskVB = (vandenberk[:,0] > 1600) & (vandenberk[:,0] < 1700) 
#            norm = np.median(vandenberk[:,1][maskVB]) / np.median(flux[mask])
#            ax.step(vandenberk[:,0],vandenberk[:,1]/norm, hold=True,label='Composite')
    
    

            if n == 'NIR':
                ax.plot(0,0,color=color, label= name)
        if errorbars == True:
            #Plotting errorbars
            ax.fill_between(wl,
                fitsfile[0].data - fitsfile[1].data, 
                fitsfile[0].data + fitsfile[1].data,
                alpha=0.4, edgecolor='blue', facecolor='blue')
        
           
            ax.fill_between(wl,
                fitsfile[0].data/trans[0].data - fitsfile[1].data/trans[0].data, 
                fitsfile[0].data/trans[0].data + fitsfile[1].data/trans[0].data,
                alpha=0.4, edgecolor='red', facecolor='red')
        
        
            ax.fill_between(wl,
                fitsfile[0].data/transdir[0].data - fitsfile[1].data/transdir[0].data, 
                fitsfile[0].data/transdir[0].data + fitsfile[1].data/transdir[0].data,
                alpha=0.4, edgecolor='green', facecolor='green')
            
    #        ax.fill_between(wl,
    #            fitsfile[0].data/transmodel - fitsfile[1].data/transmodel, 
    #            fitsfile[0].data/transmodel + fitsfile[1].data/transmodel,
    #            alpha=0.4, edgecolor='seagreen', facecolor='seagreen')
            
    from itertools import chain
    wl = np.asarray(list(chain(wlUVB,wlVIS,wlNIR)))
    flux = np.asarray(list(chain(fluxUVB,fluxVIS,fluxNIR)))
    fluxerr = np.asarray(list(chain(fluxerrUVB,fluxerrVIS,fluxerrNIR)))
    flux_corr = np.asarray(list(chain(flux_corrUVB,flux_corrVIS,flux_corrNIR)))
    flux_correrr = np.asarray(list(chain(flux_correrrUVB,flux_correrrVIS,flux_correrrNIR)))
    
    return wl, flux , fluxerr, flux_corr, flux_correrr
    
#if __name__ == '__main__':
#    ebv=0.0201
#    z=1.98041
#    offset = 1.0
#    color = None
#    printsn=False
#    sdss='spec-0284-51662-0263.fits'
#    name='SDSS1150-0023'
#    analysis(ebv,z,sdss,name,printsn,offset,color)    