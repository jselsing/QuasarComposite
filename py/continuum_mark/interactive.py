# -*- coding: utf-8 -*-


from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
import matplotlib.pyplot as plt
import methods
from pylab import pause
from matplotlib.lines import Line2D
import numpy as np

__all__ = ["continuum_mark"]


class continuum_mark(object):
    """
    A continuum estimation suite.

    Key-bindings
    
      'left-mouse-button' 
      
              place spline point

      'enter'     

              1. sort the continuum-point-array according to the x-values
              2. fit a spline and evaluate it in the wavelength points
              3. plot the continuum


      'y'
          
              Insert evenly spaced points - every 100'th pixel

      'u'
              Insert evenly spaced points - every 50'th pixel


      'd'
          
              Delete point under mouse    

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

             
        elif event.key=='enter':
            'Sort spline points and interpolate between marked continuum points'
            self.continuum, self.leg, self.con = methods.spline_interpolation(self.pointx, self.pointy, self.wave,
                               self.wave_temp, self.flux, self.flux_temp,
                               self.ax, self.leg, self.con, endpoints = self.endpoint,
                               endpoint_order = self.endpoint_order)
            self.canvas.draw() 




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

