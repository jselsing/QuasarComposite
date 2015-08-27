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


