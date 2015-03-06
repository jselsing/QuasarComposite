# -*- coding: utf-8 -*-


from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
import matplotlib.pyplot as plt
from continuum_mark import methods
from pylab import pause
from matplotlib.lines import Line2D
import numpy as np

__all__ = ["continuum_mark"]


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



def f8(seq):
    seen = set()
    return np.array([i for i, x in enumerate(seq) if x not in seen and not seen.add(x)])


class continuum_mark(object):
    """
    A normalisation suite.

    Key-bindings
    
      'left-mouse-button' 
      
              place spline point

      'enter'     

              1. sort the continuum-point-array according to the x-values
              2. fit a spline and evaluate it in the wavelength points
              3. plot the continuum
      'n'
          
              Apply normalisation
      'w'
          
              Write to file      

      'y'
          
              Insert evenly spaced points - every 100'th pixel

      'u'
              Insert evenly spaced points - every 50'th pixel


      'd'
          
              Delete point under mouse    
              
      't'
          
              Filter points by low-order chebyshev fitting and clipping values of high sigma iteratively until continuum is found        
              
      'm'
            
              Mask values (spectrum - continuum fit) > 2 sigma

      'i'
              Run through the process y-t-enter-mask until the median value of the realisations converge
    """


    def __init__(self, wave, flux, fluxerror):
        self.showverts = True
        self.epsilon = 6  # max pixel distance to count as a vertex hit
        self.pointx = []
        self.pointy = []
        self.linewidth= 0.1
        self.linewidth_over = 1.0
        self.exclude_width = 1
        self.sigma_mask = 10
        self.lover_mask = 0
        self.tolerance = 0.25
        self.leg_order = 12
        self.division = 100
        self.spacing = 100
        self.endpoint = 'n'
        self.endpoint_order = 4




        self.wave = wave
        self.flux = flux
        self.continuum = self.flux
        self.fluxerror = fluxerror
        self.wave = self.wave[~(np.isnan(self.flux))]
        self.flux = self.flux[~(np.isnan(self.flux))]
        self.fluxerror = self.fluxerror[~(np.isnan(self.flux))]
        self.ind_err = np.arange(len(self.wave))
        self.wave_temp = self.wave[self.ind_err]
        self.flux_temp = self.flux[self.ind_err]
        self.fluxerror_temp = self.fluxerror[self.ind_err]
        self.fig = plt.figure()
        canvas = self.fig.canvas
        self.ax = self.fig.add_subplot(111)
        self.line = Line2D(self.wave,self.flux, color='black', drawstyle='steps-mid',lw=self.linewidth/2,label='spectrum',zorder = 1)
        self.ax.add_line(self.line)
        self.point = None
        self.con = None
        self.leg = None
        self.llr = None
        self.diff = None
        self.error = None
        self.over = None
        
        self._ind = None # the active vert
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.left_button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.left_button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas
        
    def run(self):
        self.ax.relim()
        self.ax.autoscale_view(True,True,True)
        self.fig.canvas.manager.window.raise_()
        self.fig.canvas.draw()

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.fig.bbox)
        self.ax.draw_artist(self.line)
        if self.point != None:            
            self.ax.draw_artist(self.point)
        self.canvas.blit(self.fig.bbox)
        
    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        if self.point == None:
            ind = None
            return ind
        elif self.point != None:   
            xyt = self.ax.transData.transform(self.point.get_xydata())
            xt, yt = xyt[:, 0], xyt[:, 1]
            d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
#            dx = abs(xt-event.x)
            indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
            ind = indseq[0]
            if (d[ind]>=self.epsilon):# and dx[ind]>=self.epsilon:
                ind = None 
            return ind
            
    def mask_manual(self, mask):
        self.mask = mask     
    

     
    def left_button_press_callback(self, event):
        'whenever a left mouse button is pressed'
        if isinstance(self.pointx, np.ndarray):
            self.pointx = self.pointx.tolist()
            self.pointy = self.pointy.tolist()
        if event.inaxes==None: return
        if event.button != 1: return
        ind = self.get_ind_under_point(event)
        if ind == None:
            toolbar = plt.get_current_fig_manager().toolbar
            if event.button==1 and toolbar.mode=='':
                window = ((event.xdata-5)<=self.wave) & (self.wave<=(event.xdata+5))
                y = np.median(self.flux[window]).astype(np.float64)
                self.pointx.append(event.xdata)
                self.pointy.append(y)
                if self.point == None:
                    self.point, = self.ax.plot(self.pointx ,self.pointy ,'o' , color = 'red',ms=self.epsilon,picker=self.epsilon,label='cont_pnt')
                else:
                    self.point.remove()
                    self.canvas.draw()
                    self.point, = self.ax.plot(self.pointx ,self.pointy ,'o' , color = 'red',ms=self.epsilon,picker=self.epsilon,label='cont_pnt')
            
        self._ind = self.get_ind_under_point(event)
        self.canvas.restore_region(self.background)   
        self.canvas.draw()   
        self.canvas.blit(self.fig.bbox)

    def left_button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None
        self.canvas.draw()

    def motion_notify_callback(self, event):
        'on mouse movement'
#        if not self.showverts: return
        self._ind = self.get_ind_under_point(event)
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        toolbar = plt.get_current_fig_manager().toolbar
        if event.button==1 and toolbar.mode=='':
            x,y = event.xdata, event.ydata
            self.point.get_xydata()[self._ind] = x,y
            self.point.set_data(zip(*self.point.get_xydata()))
            self.pointx[self._ind] = x
            self.pointy[self._ind] = y
            self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.point)        
            self.canvas.blit(self.fig.bbox)

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
            
        if event.key=='d':
            'Delete point under mouse'
            self._ind = self.get_ind_under_point(event)
            if self._ind == None: return
            self.pointx = np.delete(self.pointx,[self._ind])
            # del self.pointx[self._ind]
            self.pointy = np.delete(self.pointy,[self._ind])
            # del self.pointy[self._ind]
            self.point.remove()
            self.canvas.draw()
            self.point, = self.ax.plot(self.pointx ,self.pointy ,'o' , color = 'red',ms=self.epsilon,picker=self.epsilon,label='cont_pnt')


        elif event.key=='y':
            "Insert evenly spaced points - every 100'th pixel"
            self.pointx, self.pointy, self.point = methods.insert_points(self.pointx, 
                                       self.pointy, self.wave, self.wave_temp, 
                                       self.flux, self.flux_temp, self.ax, 
                                       self.point, spacing = self.spacing / 2.0)
            self.canvas.draw()



        if event.key=='u':
            "Insert evenly spaced points - every 50'th pixel"
            self.pointx, self.pointy, self.point = methods.insert_points(self.pointx, 
                                       self.pointy, self.wave, self.wave_temp, 
                                       self.flux, self.flux_temp, self.ax, 
                                       self.point, spacing = self.spacing / 4.0)
            self.canvas.draw()
             
        elif event.key == 't':
             'Filter points by  low-order legendre fitting and clipping values of high sigma iteratively until continuum is found'
             self.pointx, self.pointy, self.point = methods.filtering(self.pointx, self.pointy, self.wave, 
                                       self.wave_temp, self.flux, self.flux_temp, 
                                       self.ax, self.point, self.leg, tolerance=self.tolerance, leg_order=self.leg_order,
                                       division=self.division)
             self.canvas.draw()
             
        elif event.key=='enter':
            'Sort spline points and interpolate between marked continuum points'
            self.continuum, self.leg, self.con = methods.spline_interpolation(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.leg, self.con, endpoints = self.endpoint,
                               endpoint_order = self.endpoint_order)
            self.canvas.draw() 
            
        elif event.key=='m':
            'Mask areas where signal is present'

            self.wave_temp, self.flux_temp, self.diff, self.error, self.over = \
                methods.mask(self.pointx, self.pointy, self.wave, 
                                       self.wave_temp, self.flux, self.fluxerror, 
                                       self.flux_temp, self.continuum, self.ax, 
                                       self.diff, self.error, self.over, 
                                       exclude_width=self.exclude_width,
                                       sigma_mask=self.sigma_mask, lower_mask_bound = self.lover_mask )
            self.canvas.draw()                            
            
        elif event.key == 'a':
            'Local linear regression'
            self.continuum, self.llr = methods.llr(self.wave, self.wave_temp, self.flux, 
                                           self.flux_temp, self.ax, self.llr)
            self.canvas.draw() 
        
        elif event.key ==  'i':
            'Iterate over points, filter, spline, mask'
            self.con_err = []
            sigma, val, i = np.std(self.flux), np.median(self.flux), 0
               
            while sigma > val / 1e3 and i < 150:
                i += 1
                if i <= 5:

                    pause(0.001)
                    self.pointx, self.pointy, self.point = methods.insert_points(self.pointx,
                               self.pointy, self.wave, self.wave_temp,
                               self.flux, self.flux_temp, self.ax,
                               self.point, spacing = self.spacing)
                    self.canvas.draw()
                    pause(0.001)
                    self.pointx, self.pointy, self.point = methods.filtering(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.point, self.leg, tolerance=self.tolerance, leg_order=self.leg_order,
                               division=self.division)
                    self.canvas.draw()
                    pause(0.001)
                    # self.pointx, self.pointy, self.point = methods.filtering(self.pointx, self.pointy, self.wave,
                    #            self.wave_temp, self.flux, self.flux_temp,
                    #            self.ax, self.point, self.leg, tolerance=self.tolerance, leg_order=self.leg_order,
                    #            division=self.division)
                    # self.canvas.draw()
                    # pause(0.001)
                    self.continuum, self.leg, self.con = methods.spline_interpolation(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.leg, self.con, endpoints = self.endpoint,
                               endpoint_order = self.endpoint_order)
                    self.canvas.draw()
                    self.con_err.append(self.continuum)

                    pause(0.001)
                    self.wave_temp, self.flux_temp, self.diff, self.error, self.over = \
                        methods.mask(self.pointx, self.pointy, self.wave,
                                           self.wave_temp, self.flux, self.fluxerror,
                                           self.flux_temp, self.continuum, self.ax,
                                           self.diff, self.error, self.over,
                                           exclude_width=self.exclude_width,
                                           sigma_mask=self.sigma_mask, lower_mask_bound = self.lover_mask  )

                else:
                    self.pointx, self.pointy, self.point = methods.insert_points(self.pointx,
                               self.pointy, self.wave, self.wave_temp,
                               self.flux, self.flux_temp, self.ax,
                               self.point, spacing = self.spacing)
                    self.pointx, self.pointy, self.point = methods.filtering(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.point, self.leg, tolerance=self.tolerance, leg_order=self.leg_order,
                               division=self.division)
                    # self.pointx, self.pointy, self.point = methods.filtering(self.pointx, self.pointy, self.wave,
                    #            self.wave_temp, self.flux, self.flux_temp,
                    #            self.ax, self.point, self.leg, tolerance=self.tolerance, leg_order=self.leg_order,
                    #            division=self.division)
                    self.continuum, self.leg, self.con = methods.spline_interpolation(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.leg, self.con, endpoints = self.endpoint,
                               endpoint_order = self.endpoint_order)
                    self.con_err.append(self.continuum)
                    self.wave_temp, self.flux_temp, self.diff, self.error, self.over = \
                        methods.mask(self.pointx, self.pointy, self.wave,
                                           self.wave_temp, self.flux, self.fluxerror,
                                           self.flux_temp, self.continuum, self.ax,
                                           self.diff, self.error, self.over,
                                           exclude_width=self.exclude_width,
                                           sigma_mask=self.sigma_mask, lower_mask_bound = self.lover_mask  )

                if i > 1:
                    sigma = abs(np.std(np.median(self.con_err, axis= 0) - np.median(self.con_err[:-1], axis= 0)))
                print i, abs(np.std(np.median(self.con_err, axis= 0) - np.median(self.con_err[:-1], axis= 0))), val/ 1e3
                self.i = i
            from gen_methods import mytotal
            self.con_err = np.array(self.con_err)
            print np.shape(self.con_err)
            from gen_methods import smooth
            self.continuum = smooth(mytotal(self.con_err, axis= 2, type='meanclip'), window_len = 1000, window = 'hanning')
            print np.shape(self.continuum)
            self.stderror = smooth(mytotal(self.con_err, axis= 2, type='stdev'), window_len = 1000, window = 'hanning')/np.sqrt(i)
            # self.continuum = np.median(self.con_err, axis= 0)
            # self.stderror = np.std(self.con_err,axis=0)/np.sqrt(i)


            self.con.remove()
            self.con = None
            continuum_mark.clear(self)

            self.canvas.draw()
            
            self.con, = self.ax.plot(self.wave,self.continuum,'r-',lw=self.linewidth_over,label='continuum',zorder = 10)

            for n in [1,2,3]:            
                self.ax.plot(self.wave,self.continuum+n*self.stderror,'r-',lw=self.linewidth_over/2.,label='continuum',zorder = 10)
                self.ax.plot(self.wave,self.continuum-n*self.stderror,'r-',lw=self.linewidth_over/2.,label='continuum',zorder = 10)

        elif event.key=='n':
            'Apply continuum_mark'
            continuum_mark.clear(self)
            self.canvas.draw()

            #self.fig.set_size_inches(14,10)
            plt.savefig(self.filename[:-5]+ "_norm.eps")

            self.flux /= self.continuum
            self.fluxerror /= self.continuum
            if hasattr(self, 'stderror'):
                self.stderror /= self.continuum
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_ylim((-0.5,1.5))
            y1 = np.ones(np.shape(self.wave))
            self.line, = self.ax.plot(self.wave,self.flux, color='black', drawstyle='steps-mid',lw=self.linewidth,label='normalised spectrum',zorder = 1)
            self.line1, = self.ax.plot(self.wave,y1, color='red', drawstyle='steps-mid',lw=self.linewidth_over,label='1',zorder = 10,alpha=1.0)
            
        elif event.key=='w':
            'Write to file'
            print 'Writing to file '+self.filename[:-5]+'_norm.fits'
            self.fitsfile[0].data = self.flux
            self.fitsfile[1].data = self.fluxerror
            self.fitsfile.writeto(self.filename[:-5]+'_norm.fits', clobber = True)
            if hasattr(self, 'stderror'):
                from astropy.io import fits
                fits.append(self.filename[:-5]+'_norm.fits', self.stderror, self.fitsfile[1].header)


        self.ax.relim()
        self.ax.autoscale_view(True,True,True)
        self.canvas.draw()

    def clear(self):
        if self.point != None:
            self.point.remove()
        self.point = None

        #if self.con != None:
        #    self.con.remove()
        #self.con = None

        if self.leg != None:
            self.leg.remove()
        self.leg = None

        if self.diff != None:
            self.diff.remove()
        self.diff = None

        if self.error != None:
            self.error.remove()
        self.error = None

        if self.over != None:
            self.over.remove()
        self.over = None

        self.pointx = []
        self.pointy = []
        self.wave_temp = []
        self.flux_temp = []

