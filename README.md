###QuasarComposite
===============
A repository to contain the the X-shooter Quasar Composite.
All code used to generate the figures and data handling are in the py-directory. The data used here is public and available under ESO program ID 090.A-0147(A). The reduction setup used is the ReflexWorkflow.kar, which uses the ESO X-shooter pipeline version 2.5.2. The reduced spectra are made available in the file py/data/regularized.dat, where the spectra are moved to rest, resampled onto a common grid and corrected for galactic extinxtion. The first column in regularized.dat is the corresponding wavelength array.

The composite itself is the file Selsing2015.dat. Feel free to download and use, but please cite as Selsing et al. 2016 A&A 585, A87. http://dx.doi.org/10.1051/0004-6361/201527096

*OBS* The composite contains zeroes where all constituent pixels are masked. This is affecting the longest wavelenghts where few spectra contribute and the sky-residuals are non-negligible.

The arxiv version is available at http://arxiv.org/abs/1510.08058v1
