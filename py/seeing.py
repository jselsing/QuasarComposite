

from __future__ import division, print_function


__author__ = 'jselsing'












if __name__ == '__main__':
    from matplotlib import rc_file
    rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
    from astropy.io import fits
    import glob
    import matplotlib.pyplot as pl
    import numpy as np

    #Files
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*SDSS*/')
    # print(sdssobjects)
    obj_list =   [ 'SDSS0820+1306', 'SDSS1150-0023', 'SDSS1219-0100', 'SDSS1236-0331' , 'SDSS1354-0013',
               'SDSS1431+0535', 'SDSS1437-0147']
    # print(sdssobjects[0][-14:-1])
    sdssobjects = [i for i in sdssobjects if i[-14:-1] in obj_list]


    arms = ['UVB', 'VIS', 'NIR']
    spat_fwhmUVB = []
    spat_fwhm_errUVB = []
    spat_fwhmVIS = []
    spat_fwhm_errVIS = []
    spat_fwhmNIR = []
    spat_fwhm_errNIR = []
    for i in sdssobjects:
        transmissionfiles = glob.glob(i+'*tran*')
        obs = glob.glob(i+'*OBJECT*/*2D*')

        obj_name = i[-14:-1]

        wl_out = []
        flux_out = []
        #test = []
        err_out = []
        print(obj_name)
        #arms = ['UVB']
        for n in arms:
            print('In arm: '+n)
            obser = [k for k in obs if n in k]

            tran = [l for l in transmissionfiles if n in l]

            ob = fits.open(obser[0])
            print(ob[0].header['CDELT1'])
            wl = 10.*(np.arange((np.shape(ob[0].data)[1]))*ob[0].header['CDELT1']+ob[0].header['CRVAL1'])

            if n == 'UVB':
                mask = (wl > 3975) & (wl < 4025)
                mu = 49
                conv = 0.16
            if n == 'VIS':
                mask = (wl > 7800) & (wl < 7850)
                mu = 49
                conv = 0.16
            if n == 'NIR':
                mask = (wl > 12975) & (wl < 13025)
                mu = 39
                conv = 0.21

            dat = np.sum(ob[1].data[:,mask], axis=1)
            if n == 'UVB':
                dat[0:33] = dat[33]
                dat[64:] = dat[64]
            if n == 'VIS':
                dat[0:31] = dat[31]
                dat[64:] = dat[64]
            if n == 'NIR':
                dat[0:29] = dat[29]
                dat[48:] = dat[48]





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
            p0V = [max(dat), mu, 5, 0, min(dat)]
            best_vals, cov = op.curve_fit(model, x_val, dat, p0 = p0V)
            # if n == "NIR":
            #     pl.plot(x_val, dat)
            #     pl.plot(x_val_pl, model(x_val_pl,*best_vals))
            #     pl.show()

            print("Seeing FWHM (pixels): ")
            fwhm = 0.5346 * best_vals[3] + np.sqrt(0.2166 * best_vals[3]**2 + best_vals[2]**2)
            dfdl = 0.5346 - 0.5 * ((0.2166 * best_vals[3]**2 + best_vals[2]**2) ** (-3/2)) *(2 * 0.2166 * best_vals[3])
            dfdg = - 0.5 * ((0.2166 * best_vals[3]**2 + best_vals[2]**2) ** (-3/2)) *(2 * best_vals[2])
            fwhm_err = np.sqrt(( dfdl**2) * (cov[3,3]**2) + ( dfdg**2) * (cov[2,2]**2))
            print(fwhm)
            print("Seeing FWHM (arcsec): ")


            print(fwhm * conv )
            print("Queried seeing (arcsec): ")
            print(ob[0].header["HIERARCH ESO TEL AMBI FWHM START"])
            print(ob[0].header["HIERARCH ESO TEL AMBI FWHM END"])
            print(ob[0].header["HIERARCH ESO TEL IA FWHMLINOBS"])
            print()
            if n == 'UVB':
                spat_fwhmUVB.append(fwhm * conv )
                spat_fwhm_errUVB.append(fwhm_err*conv)
            if n == 'VIS':
                spat_fwhmVIS.append(fwhm * conv )
                spat_fwhm_errVIS.append(fwhm_err*conv)
            if n == 'NIR':
                spat_fwhmNIR.append(fwhm * conv )
                spat_fwhm_errNIR.append(fwhm_err*conv)

    # #spat_fwhm = np.array([5.03741617638 , 4.5151829206, 4.16988640786, 4.12277668042, 4.39278362164, 4.3089644867, 5.13016814353])
    spec_fwhmVIS = 0.2 *  np.array([3.434 , 3.238, 2.663, 3.292, 3.115, 2.659, 3.322])

    spat_fwhmUVB = np.array(spat_fwhmUVB)
    spat_fwhmVIS = np.array(spat_fwhmVIS)
    spat_fwhmNIR = np.array(spat_fwhmNIR)

    spec_fwhmUVB = 4000 / ((spat_fwhmUVB / spat_fwhmVIS) * spec_fwhmVIS)
    print(np.mean(spec_fwhmUVB))
    print(np.mean(7825 / spec_fwhmVIS))
    spec_fwhmNIR = 13000 / ((spat_fwhmNIR / spat_fwhmVIS) * spec_fwhmVIS * 3)
    print(np.mean(spec_fwhmNIR))
    print(spec_fwhmNIR)
    spec_fwhmNIR = (spat_fwhmNIR / spat_fwhmVIS) * spec_fwhmVIS * 4.5
    print(spec_fwhmNIR)

    spec_fwhm = spec_fwhmVIS
    spat_fwhm = spat_fwhmVIS
    spat_fwhm_err = spat_fwhm_errVIS
    spec_fwhm_err = 0.01 * spec_fwhm

    R =  spec_fwhm
    pl.scatter(np.array(spat_fwhm), R, s= 4)

    # x = [124.46, 8.20, 52.55, 4.33]
    # y = [124.46, 50.2, 78.3, 778.8]
    # xerr = [54.2, 0.1, 2.41, 1.78]
    # yerr = [22.55, 0.37, 3.77, 0.14]
    #
    descrip = ['SDSS0820+1306', 'SDSS1150-0023', 'SDSS1219-0100', 'SDSS1236-0331', 'SDSS1354-0013', 'SDSS1431+0535', 'SDSS1437-0147']


    pl.errorbar(spat_fwhm, R, xerr=spat_fwhm_err, yerr=spec_fwhm_err, capsize=1, ls='none', color='black',
                elinewidth=1)
    s = [3 - 3 * n for n in range(len(spat_fwhm))]
    for xpos, ypos, name, yoff in zip(spat_fwhm, R, descrip, s):
        pl.annotate(name, xy =[xpos, ypos], xytext=[xpos - 25, ypos + yoff], va='bottom',
                    textcoords='offset points', fontsize=10)
    pl.title("Spatial - Spectral Seeing ")
    pl.xlabel("Spatial FWHM [arcsec]")
    pl.ylabel(r"Spectral FWHM [\AA]")
    pl.ylim((0.5, 0.8))
    pl.xlim((0.6, 0.9))
    pl.tight_layout()
    pl.savefig("../documents/figs/Seeing.pdf", dpi= 150)
    pl.show()