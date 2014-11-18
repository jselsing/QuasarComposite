

from __future__ import division, print_function


__author__ = 'jselsing'




from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')







if __name__ == '__main__':
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np

    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/')
    arms = ['UVB', 'VIS', 'NIR']
    for i in sdssobjects:
        transmissionfiles = glob.glob(i+'*tran*')
        obs = glob.glob(i+'*OBJECT*/*2D*')
        obj_name = i[-14:-1]

        wl_out = []
        flux_out = []
        #test = []
        err_out = []
        print(obj_name)
        arms = ['VIS']
        for n in arms:
#            print 'In arm: '+n
            obser = [k for k in obs if n in k]

            tran = [l for l in transmissionfiles if n in l]

            ob = fits.open(obser[0])
            print(np.shape(ob[1].data))
            dat = np.sum(ob[1].data[:,1000:-1000], axis=1)
            dat[0:31] = dat[31]
            dat[64:] = dat[64]

            def gauss(x, amp, mu, sigma):
                y = amp * np.e ** - ((( x - mu )**2) / (2*(sigma**2)))
                return y





            from scipy.special import wofz
            def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
                """
                voigt profile

                V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
                z = (x+i*gam)/(sig*sqrt(2))

                Converted from
                http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
                """
                tmp = 1.0/wofz(np.zeros((len(xarr)))+1j*np.sqrt(np.log(2.0))*Lfwhm).real
                tmp = tmp*amp*wofz(2*np.sqrt(np.log(2.0))*(xarr-xcen)/Gfwhm+1j*np.sqrt(np.log(2.0))*Lfwhm).real
                return tmp

            def model(x, amp, mu, gsigma, lsigma, b):
                y = b + voigt(x, amp, mu, gsigma, lsigma)
                return y



            import scipy.optimize as op
            x_val = np.arange(len(dat))
            x_val_pl = np.arange(0, len(dat), 0.1)
            p0G = [max(dat), 47, 5, min(dat)]
            p0V = [max(dat), 47, 5, 5, min(dat)]
            best_vals, cov = op.curve_fit(model, x_val, dat, p0 = p0V)
            #print(best_vals)
            #pl.plot(x_val, dat)
            #pl.plot(x_val_pl, model(x_val_pl,*best_vals))
            #pl.show()
            print("Seeing FWHM (pixels): ")
            fwhm = 0.5346 * best_vals[3] + np.sqrt(0.2166 * best_vals[3]**2 + best_vals[2]**2)

            print(fwhm)
            print("Seeing FWHM (arcsec): ")
            print(fwhm * 0.16)
            print("Queried seeing (arcsec): ")
            print(ob[0].header["HIERARCH ESO TEL AMBI FWHM START"])
            print(ob[0].header["HIERARCH ESO TEL IA FWHMLINOBS"])
            print()