# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:25:00 2014

@author: jselsing
"""


#================================================================================
#		Packages to use
#================================================================================



from astropy.io import fits
import glob
import Analyser
import numpy as np
#from matplotlib import rc_file
#rc_file('/Users/jselsing/Pythonlibs/plotting/qso.rc')
import matplotlib.pyplot as plt
from scipy import interpolate
from gen_methods import smooth
from scipy import stats

#================================================================================
#		Functions to be used
#================================================================================



def binning1d(array,bin):
    
    """
    Used to bin low S/N 1D  response data from xshooter.
    Calculates the biweighted mean (a la done in Kriek'10). 
    Returns binned 1dimage
    """
#    ;--------------
    s=len((array))
    outsize=s/bin
    res = np.zeros((outsize))
    for i in np.arange(0,s-(bin+1),bin):
             res[((i+bin)/bin-1)] = np.sum(array[i:i+bin])/bin
    return res



#================================================================================
#		Program
#================================================================================




#------------------------- Opsætning -------------------------
# Run Analyser on objects:
ax = None
plot = True
if plot == True:
    fig, ax = plt.subplots()
sdss1431 = Analyser.analysis(0.0274,2.096,'spec-1828-53504-0300.fits','SDSS1431+0535',plot=plot,color='black',offset = 1e+1,ax=ax)
sdss1150 = Analyser.analysis(0.0201,1.98041,'spec-0284-51662-0263.fits','SDSS1150-0023',plot=plot,color='black',offset = 1e+2,ax=ax)
sdss1236 = Analyser.analysis(0.0268,1.82389,'spec-0335-52000-0212.fits','SDSS1236-0331',plot=plot,color='black',offset = 1e+3,ax=ax)
sdss1219 = Analyser.analysis(0.0279,1.57697,'spec-0288-52000-0131.fits','SDSS1219-0100',plot=plot,color='black',offset = 1e+4,ax=ax)
sdss1354 = Analyser.analysis(0.0330,1.51237,'spec-0301-51641-0319.fits','SDSS1354-0013',plot=plot,color='black',offset = 1e+5,ax=ax)
sdss1437 = Analyser.analysis(0.0378,1.30914,'spec-0919-52409-0520.fits','SDSS1437-0147',plot=plot,color='black',offset = 1e+6,ax=ax)
sdss0820 = Analyser.analysis(0.0253,1.1257,'spec-2422-54096-0494.fits','SDSS0820+1306',plot=plot,color='black',offset = 1e+7,ax=ax)

    
    
    
if plot == True:
#    plt.yscale('log')
    plt.xlim(983,11660)
    plt.ylim(1e-18,1e-6)
    plt.yscale('log')
    plt.xlabel('Wavelength / [Angstrom]')
    plt.ylabel(r'flux / $\left[ \frac{ergs}{cm^2 \cdot sec \cdot Ang} \right]$')
    plt.legend(loc=1,prop={'size':8})
    #fig.size((8.0,5.0))
    fig.tight_layout(pad=0.1)  # Make the figure use all available whitespace
#    fig.savefig('/Users/jselsing/Documents/PhD/Conferences/XSGRB workshop - Granada/QSOcomp10.png',dpi=500)
    plt.show(block=True)
# Make returned arrays of analysed objects:
wl = [sdss0820[0],sdss1150[0],sdss1219[0],sdss1236[0],sdss1354[0],sdss1431[0],sdss1437[0]]
fext = 3
flux = [sdss0820[fext],sdss1150[fext],sdss1219[fext],sdss1236[fext],sdss1354[fext],sdss1431[fext],sdss1437[fext]]
feext = 4
fluxsig = [sdss0820[feext],sdss1150[feext],sdss1219[feext],sdss1236[feext],sdss1354[feext],sdss1431[feext],sdss1437[feext]]

# Interpolate to a common wavelength:
short = []
tall = []
for i in wl:
    short.append(min(i))
    tall.append(max(i))
short = min(short)
tall = max(tall)
step = 0.06 #CDELT
wl_new = (np.arange(((tall-short)/step))*step+short)
for n in range(np.shape(wl)[0]):    
    f = interpolate.interp1d(wl[n],flux[n],kind='linear',bounds_error = False, fill_value=0.)
    g = interpolate.interp1d(wl[n],fluxsig[n],kind='linear',bounds_error = False, fill_value=1.)
    flux[n] = f(wl_new)/np.median(f(wl_new))
    fluxsig[n] = g(wl_new)/np.median(f(wl_new))




#------------------------- Combination -------------------------
# Weighted average:
weight = 1./(np.array(fluxsig)**2)
weightarimeticmean = np.average(flux, axis = 0, weights = weight)
errofweightmean = np.sqrt(1./np.sum(np.array(fluxsig)**-2.,axis=0))
# Geometric mean:
geometricmean = stats.gmean(flux, axis = 0)
#Median:
median = np.median(flux, axis = 0)


mean = weightarimeticmean


# Calculate real S/N
sn = []
wlsn = []
step = 17
for n in range(0,len(mean),step):
#    sn.append(np.mean((mean)[n:n+step])/np.std((mean)[n:n+step]))
    sn.append(np.mean((mean/errofweightmean)[n:n+step]))
    wlsn.append(np.mean((wl_new)[n:n+step]))
    
# Print estimated S/N
print 'Estimated S/N in: (wl > 19950) & (wl < 20250)'
mask = (wl_new > 1460) & (wl_new < 1470)
print np.mean((mean)[mask])/np.std((mean)[mask])
print 'Real S/N in: (wl > 19950) & (wl < 20250)'
mask = (wl_new > 1460) & (wl_new < 1470)
print np.mean((mean)[mask])/np.mean((errofweightmean)[mask])



#------------------------- Auxilliary products -------------------------
# Johans Composite:
vandenberk = np.loadtxt('compom.data')
mask = (wl_new > 1600) & (wl_new < 1700) 
maskVB = (vandenberk[:,0] > 1600) & (vandenberk[:,0] < 1700) 
norm = np.median(vandenberk[:,1][maskVB]) / np.median(mean[mask])

# Vanden Berk:
vandenberk = np.loadtxt('composite.txt')
mask = (wl_new > 1600) & (wl_new < 1700) 
maskVB = (vandenberk[:,0] > 1600) & (vandenberk[:,0] < 1700) 
norm = np.median(vandenberk[:,1][maskVB]) / np.median(mean[mask])


#------------------------- Plotting -------------------------

# Different combination methods:
bin = 1
fig, ax = plt.subplots()

ax.plot(binning1d(wl_new,bin)[:-10],binning1d(mean,bin)[:-10],drawstyle='steps-mid',color='black',label='Weighted Arithmeric Mean Composite',linewidth=0.3)


#ax.fill_between(wl_new,
#    mean - errofweightmean, 
#    mean + errofweightmean,
#    alpha=0.4, edgecolor='blue', facecolor='blue') 

#ax.plot(binning1d(wl_new,bin),binning1d(geometricmean,bin),drawstyle='steps-mid',color='red',label='Geometric Mean Composite',linewidth=1)

#ax.plot(binning1d(wl_new,bin),binning1d(median,bin),drawstyle='steps-mid',color='green',label='Median Composite',linewidth=1)

# S/N 
#ax.plot(wlsn,smooth(np.array(sn),window_len=1, window = 'hanning'))



#step = (tall-short)/10.0
#for n in range(1,3):
#    for i in range(1,6):
#        print i+(5*(n-1))
#
##        ax.subplot(5,2,i+(5*(n-1)))
#        ax.add_subplot(5,2,i+(5*(n-1)))
#        ax.plot(wl_new,mean,drawstyle='steps-mid',color='black',label='New Composite Direct Correction',linewidth=1)
#        ax.fill_between(wl_new,
#            mean - errofweightmean, 
#            mean + errofweightmean,
#            alpha=0.4, edgecolor='blue', facecolor='blue')        
#        
#        # Overplot Johans composite
#        ax.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='CompoM Composite')
##        
#        # OVerplot Vanden berk
#        ax.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Vanden Berk Composite')
#        
#        ##Overplot identified lines from composite
#        composite = np.genfromtxt('Vandenberkcompositelinelist.txt', dtype=None)
#        for k in range(len(composite)):
#            ax.axvline(x=composite[k,][1],color='green',linestyle='dashed', hold=True)
#            mask2 = (wl_new > composite[k,][1] - 5) & (wl_new < composite[k,][1] + 5)
#            val = np.median(mean[mask2])
#            ax.annotate(composite[n,][0],xy=(composite[k,][1],val*1.2),fontsize='x-small')
#        
#        #Set ranges 
#        mask = (wl_new > short + step * (i+(5*(n-1)) - 1)) & (wl_new < short + step * (i+(5*(n-1))))
#        maxv = max(mean[mask])
#        if i+(5*(n-1)) == 1:
#            maxv = 0.05*maxv
#        minv = min(mean[mask])
#        ax.set_xlim(short + step * (i+(5*(n-1)) - 1),short + step * (i+(5*(n-1))))
#        ax.set_ylim(0.9*minv,1.1*maxv)
#        plt.yscale('log')
        # Overplot Johans composite
ax.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='CompoM Composite')
#        
        # OVerplot Vanden berk
ax.plot(vandenberk[:,0],vandenberk[:,1]/norm,drawstyle='steps-mid',label='Vanden Berk Composite')
        
        ##Overplot identified lines from composite
composite = np.genfromtxt('Vandenberkcompositelinelist.txt', dtype=None)
for k in range(len(composite)):
    ax.axvline(x=composite[k,][1],color='green',linestyle='dashed')
    mask2 = (wl_new > composite[k,][1] - 5) & (wl_new < composite[k,][1] + 5)



#plotting options:
plt.title('QSOs')
plt.xlim(983,11660)
#plt.ylim(1e-18,1e-6)
plt.yscale('log')
plt.xlabel('Wavelength / [Angstrom]')
plt.ylabel(r'normalized flux')

plt.legend(loc=1,prop={'size':8})
fig.tight_layout(pad=0.1)  # Make the figure use all available whitespace
fig.savefig('/Users/jselsing/Documents/PhD/Conferences/XSGRB workshop - Granada/QSOcomp13.png',dpi=500)
plt.show()


hej = fits.open('/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/Products_final/SDSS0820+1306/FLUX_SLIT_FLUX_MERGE1D_UVB.fits')
hej[0].header['CDELT1'] = 0.06
hej[0].header['CRVAL1'] = short

fits.writeto('comp.fits',mean, hej[0].header, clobber=True)  


