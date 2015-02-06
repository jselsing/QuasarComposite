# -*- coding: utf-8 -*-
"""
Short script to plot the accuracy of the optimal template for the telluric correction
"""


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
    import numpy as np
    import glob
    import plotting.plotting as plot
    import matplotlib.pyplot as pl
    
    root_dir = '/Users/jselsing/Work/X-Shooter/CompositeRedQuasar/processed_data/'
    sdssobjects = glob.glob(root_dir+'*/')
    print(sdssobjects)
    arms = ['UVB', 'VIS', 'NIR']     
    for i in sdssobjects:
        obs = glob.glob(i+'*QC*')
        obs_check = glob.glob((i + '*uncorrected*'))

        obj_name = i[-14:-1]

#        fig = None

        if obs != []:
            for n in arms:
#                print 'In arm: '+n
                
                ob = [k for k in obs if n in k]


                dat = np.genfromtxt(str(ob[0]), dtype = np.float64)
                wl = dat[:,0]
                tell = dat[:,1]
                fit = dat[:,2]
                corr = fit/tell

                dat_check = np.genfromtxt(str(obs_check[0]), dtype = np.float64)
                wl_check = dat_check[:,0]
                flux_check = dat_check[:,1] / np.median(dat_check[:,1])

                fig = plot.plot_data(wl, corr, lw=0.2)
                fig = plot.plot_data(wl_check, flux_check, lw=0.6, color = "red", fig=fig, title = 'Hip040217', ylabel = 'Normalised Flux')
                # fig = plot.plot_data(wl,fit, lw=0.6, color = "red", fig=fig, title = 'Hip040217', ylabel = 'Normalised Flux')
                if n == 'VIS':
#                    print fig.axes                    
#                    .set_xlim((6200,10130))
                    fig.savefig('tell_corr_QC_'+n+'.eps', format='eps')
            pl.show(block=True)
            print fig.axes
                
