#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SYNOPSIS

    TODO helloworld [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    TODO: Name <name@example.org>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

VERSION

    $
"""

from __future__ import division, print_function

__all__ = ["voigt", "gauss", "continuum_fit", "wavelength_conversion", "common_wavelength"]
__version__ = "0.0.1"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2015 Jonatan Selsing"

import numpy as np

def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
    """
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Converted from
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    """
    from scipy.special import wofz
    tmp = 1.0/wofz(np.zeros((len(xarr)))+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    tmp = tmp*amp*wofz(2*np.sqrt(np.log(2.0))*(xarr-xcen)/Gfwhm+1j*np.sqrt(np.log(2.0))*Lfwhm).real
    return tmp

def gauss(xarr, amp, mu, sig2):
    """
    Generic gaussian function
    :param xarr: x-array
    :param amp: amplitude
    :param mu: x-center
    :param sig2: variance
    :return: gaussian function
    """
    tmp = amp * np.exp(-0.5 * (xarr - mu) ** 2 / sig2)
    return tmp


def linear(xarr, a, b):
    """
    :param xarr: x-array
    :param a: slope
    :param b: intercept
    :return: linear function
    """
    tmp = a * x + b
    return tmp

def power(xarr, a, b):
    """
    :param xarr: x-array
    :param a: power
    :param b: amplitude
    :return: powerlaw function
    """
    tmp = b * (xarr ** a)
    return tmp

def continuum_fit(wl, flux, fluxerr, edge_mask_len=50):
    """
    Small function to estimate continuum. Takes spectrum and uses only the edges, as specified by edge_mask_len, fits
    a polynomial and returns the fit.
    :param wl: Wavelenght array in which to estimate continuum
    :param flux: Corresponding flux array
    :param fluxerr: Corresponding array containing flux errors
    :param edge_mask_len: Size of edges to use to interpolate continuum
    :return: Interpolated continuum
    """
    from numpy.polynomial import chebyshev
    wl_chebfit = np.hstack((wl[:edge_mask_len],wl[-edge_mask_len:]))
    mask_cont = np.array([i for i,n in enumerate(wl) if n in wl_chebfit])
    chebfit = chebyshev.chebfit(wl_chebfit, flux[mask_cont], deg = 1, w=1/fluxerr[mask_cont]**2)
    chebfitval = chebyshev.chebval(wl, chebfit)
    return chebfitval, chebfit

def wavelength_conversion(input_wavelenght, conversion = 'vacuum_to_air'):
    """Convert wavelenghts between air and vacuum.
    Coversion from vacuum to air, taken from: https://www.sdss3.org/dr9/spectro/spectro_basics.php"""
    if conversion == 'vacuum_to_air':
        output_wavelenght = input_wavelenght / (1.0 + 2.735182E-4 + 131.4182 / input_wavelenght**2 + 2.76249E8 / input_wavelenght**4)

    return output_wavelenght

#Nearest naeightbour positive definite
import numpy as np,numpy.linalg

def _getAplus(A):
    eigval, eigvec = np.linalg.eig(A)
    Q = np.matrix(eigvec)
    xdiag = np.matrix(np.diag(np.maximum(eigval, 0)))
    return Q*xdiag*Q.T

def _getPs(A, W=None):
    W05 = np.matrix(W**.5)
    return  W05.I * _getAplus(W05 * A * W05) * W05.I

def _getPu(A, W=None):
    Aret = np.array(A.copy())
    Aret[W > 0] = np.array(W)[W > 0]
    return np.matrix(Aret)

def nearPD(A, nit=10):
    n = A.shape[0]
    W = np.identity(n)
# W is the matrix used for the norm (assumed to be Identity matrix here)
# the algorithm should work for any diagonal W
    deltaS = 0
    Yk = A.copy()
    for k in range(nit):
        Rk = Yk - deltaS
        Xk = _getPs(Rk, W=W)
        deltaS = Xk - Rk
        Yk = _getPu(Xk, W=W)
    return Yk

def nearPSD(A,epsilon=0):
   n = A.shape[0]
   eigval, eigvec = np.linalg.eig(A)
   val = np.matrix(np.maximum(eigval,epsilon))
   vec = np.matrix(eigvec)
   T = 1/(np.multiply(vec,vec) * val.T)
   T = np.matrix(np.sqrt(np.diag(np.array(T).reshape((n)) )))
   B = T * vec * np.diag(np.array(np.sqrt(val)).reshape((n)))
   out = B*B.T
   return(out)

def ConfInt(x,func,coeff,cov,p,nsample=5000):
  """
  Monte Carlo generation of confidence interval.
  Given the variables 'coeff' and their covariance matrix 'cov', return
  for each bin in the x-axis 'x' the confidence intervals of the function
  'func' defined by the 'p'th percentiles (single number or list of numbers),
  calculated from 'n' random generations.

  TODO: Currently, p must be a list of two numbers. Should be any number.
  """
  #funcdict = {'power':  power,
  #            'linear': linear,
  #            'gauss': gauss,
  #            'voigt': voigt}
  #f    = funcdict[func]
  f = func
  ndim = len(coeff)
  nx   = len(x)
  rlz  = np.ndarray(shape=(nsample,nx))   #Array for the nsample realizations
  import warnings
  warnings.simplefilter("error", RuntimeWarning)
  try:
      args = np.random.multivariate_normal(coeff,cov,nsample).T
  except RuntimeWarning:
      print('Warning: Covariance matrix is not positive-definite. Off-diagonal elements set to zero')
      cov = np.array(np.diagflat(np.diag(cov)))
      args = np.random.multivariate_normal(coeff,100*cov,nsample).T

  np.random.seed(1234325)
  for i in range(nsample):
    perc     = args[:,i]
    rlz[i,:] = f(x,*perc)                 #i'th realization

  lo = np.empty_like(x)
  hi = np.empty_like(x)
  for i in range(nx):
    lo[i],hi[i] = np.percentile(rlz[:,i],p)

  return lo,hi
#------------------------------------------------------------------------------



def common_wavelength(wlarr_old, wlarr_new, fluxarr_old, fill_value = 0.):
    from scipy import interpolate
    f = interpolate.interp1d(wlarr_old, fluxarr_old, kind='linear', bounds_error = False, fill_value=fill_value)
    fluxarr_new = f(wlarr_new)
    return fluxarr_new


import matplotlib
from math import sqrt
SPINE_COLOR = 'gray'

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + str(fig_height) +
              "so will reduce to " + str(MAX_HEIGHT_INCHES) + "inches.")
        fig_height = MAX_HEIGHT_INCHES
    size = 8
    params = {'backend': 'Qt4Agg',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': size, # fontsize for x and y labels (was 10)
              'axes.titlesize': size,
              'font.size': size, # was 10
              'legend.fontsize': size, # was 10
              'xtick.labelsize': size,
              'ytick.labelsize': size,
              # 'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)


def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(True)

    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax


# from gen_methods import medfilt, smooth
def hist(rawData,xRange,nBins=10,mode='lin'):

    """histogram using linear binning of supplied data

    Input:
    rawData	-- list containing data to be binned
    xRange  -- lower(incl)/upper(excl) boundary for numerical values
    nBin    -- desired number of bins (default =10)
    mode	-- binning type (possible choices: lin, log)

    Returns: (nothing)
    """
    from math   import sqrt,floor,log,exp
    h = [0]*nBins
    xMin=float(xRange[0])
    xMax=float(xRange[1])

    if mode == 'lin':
        dx = (xMax-xMin)/nBins
        def binId(val):   return int(floor((val-xMin)/dx))
        def bdry(bin):	  return xMin+bin*dx, xMin+(bin+1)*dx
        def GErr(q,n,dx): return sqrt(q*(1-q)/(N-1))/dx

    elif mode == 'log':
        dx = log(xMax/xMin)/nBins
        def binId(val):   return int(floor(log(val/xMin)/dx))
        def bdry(bin):	  return xMin*exp(bin*dx), xMin*exp((bin+1)*dx)
        def GErr(q,n,dx): return "##"

    for value in rawData:
        if 0<=binId(value)<nBins:
          h[binId(value)] += 1

    N=sum(h)
    binned = []
    for bin in range(nBins):
        hRel   = float(h[bin])/N
        low,up = bdry(bin)
        binned.append(low)
        width  = up-low
        # print(low, up, hRel/width, GErr(hRel,N,width))
    return binned




def test():
    """ Testing Docstring"""
    pass


if __name__ == '__main__':
    test()