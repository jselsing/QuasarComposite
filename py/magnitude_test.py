


import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import glob
from methods import common_wavelength
from gen_methods import medfilt
from astropy.cosmology import Planck13 as cosmo
import sys
sys.path.append('/Users/jselsing/github/pysysp/pysysp/')
import pysysp



def main():
    fil = pysysp.StarSpectrum()
    fil.setwavelength(np.arange(1000, 30000, 0.4))
    fil.setflux( 1e-15 * (fil.wavelength/10000.0)**(-1.5))
    i = pysysp.BandPass('/Users/jselsing/github/pysysp/pysysp/filters/ugriz/i.dat')


    dl = (cosmo.luminosity_distance(1.2)) * 1e5
    Mz0 = -5 * np.log10(dl.value) + fil.apmag(band=i, mag='AB')
    print('Mz0 = '+str(Mz0))


    i.wavelength /= (1 + 2)
    # fil.flux = fil.flux * (1 + 2)   

    i.update_bandpass()

    dl = (cosmo.luminosity_distance(1.2)) * 1e5
    Mz2 = -5 * np.log10(dl.value) + fil.apmag(band=i, mag='AB') #- 2.5 * np.log10(1 + 2)
    print('Mz2 = '+str(Mz2))




    # i.wavelength *= (1 + 2)

    # i.update_bandpass()

    # fil.wavelength = fil.wavelength * (1 + 2)
    # fil.flux = fil.flux / (1 + 2)
    # dl = (cosmo.luminosity_distance(1.2)) * 1e5
    # Mz2 = -5 * np.log10(dl.value) + fil.apmag(band=i, mag='AB')
    # print('Mz2 = '+str(Mz2))
    print('Mz0 - Mz2 = '+str(Mz0 - Mz2))

if __name__ == '__main__':
	main()