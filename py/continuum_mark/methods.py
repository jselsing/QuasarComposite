# -*- coding: utf-8 -*-

__all__ = ["insert_points", "llr", "spline_interpolation", "filtering", "mask"]



import numpy as np
import types


#=======================================----------------------------------------
def f8(seq):
    """Returns array with unique values"""
    seen = set()
    return np.array([i for i, x in enumerate(seq) if x not in seen and not seen.add(x)])
#=======================================----------------------------------------



#=======================================----------------------------------------
def insert_points(pointx, pointy, wave, wave_temp, flux, flux_temp, axis, point_instance, 
                  window_width = 15, spacing = 100,  pick_epsilon = 6):
    """Insert evenly spaced points"""
    if isinstance(pointx, np.ndarray):
        pointx = pointx.tolist()
        pointy = pointy.tolist()
    if not isinstance(point_instance, types.NoneType):
        pointx = []
        pointy = []
        point_instance.remove()
        point_instance, = axis.plot(pointx, pointy, 'o', color = 'red', ms = pick_epsilon, 
                           picker = pick_epsilon, label='cont_pnt')

    for i in wave_temp[30+np.random.randint(30):-30+np.random.randint(30):spacing + np.random.randint(spacing/10.0)]:
        window = np.arange(np.where(wave_temp == i)[0]-window_width,np.where(wave_temp == i)[0]+window_width)
        window2 = np.arange(np.where(wave == i)[0]-window_width,np.where(wave == i)[0]+window_width)
        if max(wave_temp[window]) - min(wave_temp[window]) == max(wave[window2]) - min(wave[window2]):
            y_val = np.median(flux_temp[window]).astype(np.float64)
            pointx.append(i)
            pointy.append(y_val)
            try:
                point_instance.remove()
            except AttributeError:
                pass
            finally:
                point_instance, = axis.plot(pointx, pointy, 'o', color = 'red', ms=pick_epsilon, 
                                   picker = pick_epsilon, label='cont_pnt')

    return np.array(pointx), np.array(pointy), point_instance
#=======================================----------------------------------------



#=======================================----------------------------------------
def llr(wave, wave_temp, flux, flux_temp, axis, llr_instance, linewidth = 2.0):
    """Local linear regression"""
    import pyqt_fit.kernel_smoothing as smooth
    
    k1 = smooth.LocalLinearKernel1D(wave_temp, flux_temp)
#    k1 = smooth.LocalPolynomialKernel1D(wave_temp, flux_temp, q= 2)
    continuum = k1(wave)
    try:
        llr_instance.remove()
    except AttributeError:
        pass
    finally:
           llr_instance, = axis.plot(wave, continuum, color = 'purple', lw = linewidth, 
                            label = 'legendre fit', zorder = 10)

  
    return continuum, llr_instance
#=======================================----------------------------------------


#=======================================----------------------------------------        
def spline_interpolation(pointx, pointy, wave, wave_temp, flux, flux_temp, axis,
                         leg_instance, con_instance, linewidth= 2.0, endpoints = 'y',
                         endpoint_order = 4):
    """Sort spline points and interpolate between marked continuum points"""
    from numpy.polynomial import chebyshev
    from scipy.interpolate import splrep,splev    
    
    #Insert endpoints
    if endpoints == 'y':
        sort_array = np.argsort(pointx)
        x = np.array(pointx)[sort_array]
        y = np.array(pointy)[sort_array]
        chebfit = chebyshev.chebfit(x, y, deg = endpoint_order)
        chebfitval = chebyshev.chebval(wave, chebfit)
        i =wave[150], wave[-150]
        window1, window2 = ((i[0]-70)<=wave) & (wave<=(i[0]+70)), ((i[1]-70)<=wave) & (wave<=(i[1]+70))
        y_val = np.median(chebfitval[window1]).astype(np.float64),np.median(chebfitval[window2]).astype(np.float64)
        pointx = np.concatenate([pointx,i])
        pointy = np.concatenate([pointy,y_val])
        ind_uni = f8(pointx)
        pointx = np.array(pointx)[ind_uni]
        pointy = np.array(pointy)[ind_uni]
        try:
            leg_instance.remove()
        except AttributeError:
            pass
        finally:
               leg_instance, = axis.plot(wave, chebfitval, color='black', lw = linewidth,
                                         label = 'legendre fit', zorder = 10)


    # print(pointx,pointy)
    #Sort numerically
    sort_array = np.argsort(pointx)
    x = np.array(pointx)[sort_array]
    y = np.array(pointy)[sort_array]


    from gen_methods import smooth
    #Interpolate
    spline = splrep(x,y, k=3)
    continuum = splev(wave,spline)
    continuum = smooth(continuum, window_len=1, window='hanning')
    
    #Plot
    try:
        con_instance.remove()
    except AttributeError:
        pass
    finally:
           con_instance, = axis.plot(wave,continuum, color='red', lw = linewidth, 
                                     label = 'continuum', zorder = 10)  
                                     
    return continuum, leg_instance, con_instance

#=======================================----------------------------------------



#=======================================----------------------------------------
def filtering(pointx, pointy, wave, wave_temp, flux, flux_temp, axis,
                         point_instance, leg_instance, linewidth= 2.0, 
                         pick_epsilon = 6, tolerance =0.05, leg_order = 1, division = 50):
    'Filter points by  low-order legendre fitting and clipping values of high sigma iteratively until continuum is found'
    # import pyqt_fit.kernel_smoothing as smooth
    from gen_methods import smooth
    #Ensure uniqueness
    ind_uni = f8(pointx)
    pointx = np.array(pointx)[ind_uni]
    pointy = np.array(pointy)[ind_uni]        
    
    #Fit with chebyshev polynomial and clip point furthest from the fit to remove points not varying smoothly
    from numpy.polynomial import chebyshev

    sigma = np.average(pointy)*tolerance
    chebfitval = 0
    while max(np.sqrt((pointy - chebfitval)**2)) >= abs(sigma):
        sort_array = np.argsort(pointx)
        x = np.array(pointx)[sort_array]               
        y = np.array(pointy)[sort_array]
        # chebfit = chebyshev.chebfit(x, y, deg = leg_order)
        # chebfitval = chebyshev.chebval(x, chebfit)
        chebfitval = smooth(x, window_len=3, window='hanning')

        ind = [i for i, j in enumerate(np.sqrt((pointy - chebfitval)**2)) if j == max(np.sqrt((pointy - chebfitval)**2))]
        pointx = np.delete(pointx,[ind[0]])
        pointy = np.delete(pointy,[ind[0]])
        chebfitval = np.delete(chebfitval,[ind[0]])

    min_sep = [min(abs(np.array(pointx)[np.nonzero(abs(np.array(pointx) - j))] - j)) for i, j in enumerate(pointx)]

    while min(min_sep) <= (max(wave) - min(wave))/(division):
        ind_min = np.where(min(min_sep) == min_sep)
        p = np.random.randint(len(ind_min[0]))
        pointx = np.delete(pointx,[(ind_min[0])[p]])
        pointy = np.delete(pointy,[(ind_min[0])[p]])
        min_sep = [min(abs(np.array(pointx)[np.nonzero(abs(np.array(pointx) - j))] - j)) for i, j in enumerate(pointx)]

    if len(pointx) >= 6:
        sort_array = np.argsort(pointx)
        x = np.array(pointx)[sort_array]
        y = np.array(pointy)[sort_array]
        # chebfit = chebyshev.chebfit(x, y, deg = leg_order)
        # chebfitval = chebyshev.chebval(x, chebfit)
        chebfitval = smooth(x, window_len=3, window='hanning')
        ind = [i for i, j in enumerate(np.sqrt((pointy - chebfitval)**2)) if j == max(np.sqrt((pointy - chebfitval)**2))]
        pointx = np.delete(pointx,[ind[0]])
        pointy = np.delete(pointy,[ind[0]])

    point_instance.remove()
    point_instance, = axis.plot(pointx ,pointy ,'o' , color = 'red',ms=pick_epsilon,picker=pick_epsilon,label='cont_pnt')

    return pointx, pointy, point_instance
#=======================================----------------------------------------


    
#=======================================----------------------------------------   
def mask(pointx, pointy, wave, wave_temp, flux, fluxerror, flux_temp, continuum, 
         axis, diff_instance, error_instance, over_instance, linewidth= 0.2, 
         exclude_width = 20, sigma_mask = 3, lower_mask_bound = 0):
    """Mask areas where signal is present"""
#        import pyqt_fit.kernel_smoothing as smooth
    from gen_methods import smooth
    difference = smooth(abs(continuum - flux), window_len=20, window='hanning')
#        k1 = smooth.SpatialAverage(self.wave, self.continuum - self.flux)
#        self.difference = k1(self.wave)
    if diff_instance != None and error_instance != None:
        diff_instance.remove()
        error_instance.remove()
    diff_instance, = axis.plot(wave, difference, color='green', lw = linewidth, 
                               label='continuum', drawstyle = 'steps-mid')
    error_instance, = axis.plot(wave, sigma_mask*fluxerror, 'blue', lw = linewidth/2.0, 
                                label = 'continuum')


    ind_err = np.array([i for i, j in enumerate(fluxerror) if difference[i] < sigma_mask*j and flux[i] >= lower_mask_bound])# and 6*fluxerror[i] < flux[i]]) #and np.average(self.difference[i-100:i+100]) <= j] )
    b = exclude_width
    ind_err = np.array([j for i, j in enumerate(ind_err[:-b]) if j + b == ind_err[i+b] and j - b == ind_err[i-b]])
    wave_temp = wave[ind_err]
    flux_temp = flux[ind_err]
    
    try:
        over_instance.remove()
    except AttributeError:
        pass
    over_instance, = axis.plot(wave_temp, flux_temp, color = 'red', lw=linewidth/4.0)
    
    return wave_temp, flux_temp, diff_instance, error_instance, over_instance
#=======================================----------------------------------------



        