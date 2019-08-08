#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import re
from numpy.ma import is_masked
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from numpy.linalg import eig, inv
from astropy.stats import sigma_clip
from . import yanny

import mangatools as mangatools_package 
package_path = os.path.dirname(os.path.realpath(mangatools_package.__file__))

def processBar(total, progress, barlength=40):
    barlength, status = barlength, ""
    progress = float(progress) / float(total)
    if progress >= 1:
        progress, status = 1, "\r\n"
    block = int(round(barlength * progress))
    text = "\r[{}] {:.0f}% {}".format("#" * block + 
                                      "-" * (barlength - block),
                                      round(progress * 100, 0), status)
    sys.stdout.write(text)
    sys.stdout.flush()

def sector_binning(r, phi, r_min=0, r_max=np.inf, phi_min=0, phi_max=360):
    """generate sector binning by elliptical coordinate
    
    Args:
        r: the radius matrix
        phi: the pich angle of the matrix
        r_min, r_max: the range of selected radius
        phi_min, phi_max: the range of selected phi in degree
    Returns:
        return the selecting ndarray with selected region with True value
    """
    return ((r < r_max) & (r >= r_min)) & ((phi < phi_max) & (phi >= phi_min))

def torus_binning(x, y, b2a, phi, Rmin=None, Rmax=None, binNum=None):
    """generate the torus_binning

    Like voronoi 2d binning, this function generate the 2d binning used for
    explore diffrent stellar type torus around the specified core

    Args:
        x: the x coordinates of given galaxy, 2d shape(1d maybe, not test)
        y: same as x, must with the same shape with x
        b2a: the rate bettween the semi-minor and semi-major axis of the ellipse
        phi: the angle which semi-major axis deviate from x coordinates
             (counter-clockwise direction)
        Rmin: the minimal R of the torus
        Rmax: the maximum R of the torus
        binNim: the flag added to the selected region

    Returns:
        binMatrix: shape like x and y, contain the same binning number for binned
                   regions
    """
    assert x.shape == y.shape, "Input matrix mis-match"
    # change to new coordinate system(counter-clockwise rotate)
    x_new = x*np.cos(phi) + y*np.sin(phi)
    y_new = y*np.cos(phi) - x*np.sin(phi)
    binMatrix = np.zeros(x.shape)
    if Rmin == None or Rmin == 0:
        a_max, b_max = Rmax, b2a*Rmax
        binMatrix[(x_new/a_max)**2 + (y_new/b_max)**2 < 1] = binNum
    else:
        a_min, b_min = Rmin, b2a*Rmin
        a_max, b_max = Rmax, b2a*Rmax
        binMatrix[((x_new/a_max)**2 + (y_new/b_max)**2 < 1) & 
               ((x_new/a_min)**2 + (y_new/b_min)**2 > 1)] = binNum
    return binMatrix

def binning(signal, binMatrix):
    """binning data with given bin-number

    Args:
        signal: the signal need to bin
        binMatrix: the binning matrix used for binning

    Return:
        binned signal, same size like signal
    """
    signal_with_bins = signal.copy() # void to destroy origin data
    for bin_num in np.unique(binMatrix):
        chose_bin = (binMatrix == bin_num)
        signal_with_bins[chose_bin] = signal[chose_bin].mean()
    return signal_with_bins

def check_in_ellipse(x, y, a, b, phi):
    # return the number of points in the ellipse
    # phi: the angle related to the x-axis direction
    x_new = x*np.cos(phi) + y*np.sin(phi)
    y_new = y*np.cos(phi) - x*np.sin(phi)
    return np.sum((x_new/a)**2 + (y_new/b)**2 < 1)

def emline_dictionary(hdu, ext):
    """get the emission line from the MaNGA DAP header
    Args:
        hdu: the header data unit
        ext: the name of the table contains multiple channel, for example 
             the emline in hdu 
    Return:
        emline dictionary, with the name of the emission line and key of it's 
        index
    """
    emlinedict = {}
    for name, value in hdu[ext].header.items():
        if name[0] == 'C':
            try:
                i = int(name[1:])-1
            except ValueError:
                continue
            emlinedict[value] = i
    return emlinedict

def bptregion(x, y, mode='N2'):
    '''
    Args:
        lines: dictionary contains all the needed for bpt
            ['Hb-4862', 'OIII-5008','Ha-6564','NII-6585']
      x: log10(NII/Ha) or log10(SII/Ha) or log10(OI/Ha)
      y: log10(OIII/Hb)
      mode: mode "N2" -> x = log10(NII/Ha)
            mode "S2" -> x = log10(SII/Ha)
            mode "O1" -> x = log10(OI/Ha)  ! not surpport yet
    Note:
      x, y should be masked array, 
      example: x = np.ma.array(x)
    '''
    # check starforming, composite or AGN
    # region = np.zeros_like(lines[0].data)
    if mode == 'N2':
        ke01 = 0.61/(x-0.47)+1.19
        ka03 = 0.61/(x-0.05)+1.3
        schawinski_line = 1.05*x+0.45
        region_AGN = np.logical_or(np.logical_and(x<0.47, y>ke01), x>0.47)
        region_composite = np.logical_and(y<ke01, y>ka03)
        region_starforming = np.logical_and(x<0.05, y<ka03)
        # depleted
        #region_seyfert = np.logical_and(x>np.log10(0.6), y>np.log10(3.))
        #region_liner = np.logical_and(region_AGN, np.logical_and(x>np.log10(0.6), y<np.log10(3.)))
        # adapted from Schawinski2007
        region_seyfert = np.logical_and(region_AGN, y>schawinski_line)
        region_liner = np.logical_and(region_AGN, y<schawinski_line)
        if is_masked(x) or is_masked(y):
            return region_AGN.filled(False), region_composite.filled(False), region_starforming.filled(False), region_seyfert.filled(False), region_liner.filled(False)
        else:
            return region_AGN, region_composite, region_starforming, region_seyfert, region_liner

    if mode == 'S2':
        ke01_line = 0.72/(x-0.32)+1.3
        seyfert_liner_line = 1.89*x+0.76
        region_seyfert = np.logical_and(np.logical_or(y>ke01_line, x>0.32), y>seyfert_liner_line)
        region_liner = np.logical_and(np.logical_or(y>ke01_line, x>0.32), y<seyfert_liner_line)
        region_starforming = np.logical_and(y<ke01_line, x<0.32)
        if is_masked(x) or is_masked(y):
            return region_seyfert.filled(False), region_liner.filled(False), region_starforming.filled(False)
        else:
            return region_seyfert, region_liner, region_starforming

def bar_statistic(data, title=None, index=None, fig=None, showImage=True):
    # data can be a list of string or number
    if index == None:
        index = 111
    if fig == None:
        fig = plt.figure()
    ax = fig.add_subplot(index)
    items, counts = np.unique(data, return_counts=True)
    total = np.sum(counts)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = np.arange(0, len(items))
    ax.bar(index, counts, 0.5, label=items)
    for i, j in enumerate(counts):
        ax.text(i-0.1, j+1, '{:.2f}%'.format(j/total*100))
    ax.set_xticks(index)
    ax.set_xticklabels(items)
    if title:
        ax.set_title(title)
    if showImage:
        plt.show()

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    ## the original algrithm
    # E, V =  eig(np.dot(inv(S), C))
    # n = np.argmax(np.abs(E))
    # a = V[:,n]
    tmpA = S[0:3,0:3]
    tmpB = S[0:3,3:6]
    tmpC = S[3:,3:]
    tmpD = C[0:3,0:3]
    tmpE = np.dot(inv(tmpC), tmpB.conj().T)
    tmpF = np.dot(tmpB,tmpE)
    eval_x, evec_x = eig(np.dot(inv(tmpD), (tmpA - tmpF)))

    # Find the negative (as det(tmpD) < 0) eigenvalue
    I = np.argmax(eval_x)
  
    #Extract eigenvector corresponding to negative eigenvalue
    a = evec_x[:,I]
  
    # Recover the bottom half...
    evec_y = np.dot(-1*tmpE, a)
    a = np.concatenate((a,evec_y))  

    return a

def get_ellipse_metrics( a ):
    thtarad = .5 * np.arctan2(a[1], a[0]-a[2])
    cost = np.cos(thtarad)
    sint = np.sin(thtarad)
    sinsq = sint*sint
    cossq = cost*cost
    cossin = sint*cost
    
    Ao = a[5]
    Au = a[3]*cost + a[4]*sint
    Av = -1*a[3]*sint + a[4]*cost
    Auu = a[0]*cossq + a[2]*sinsq + a[1]*cossin
    Avv = a[0]*sinsq + a[2]*cossq - a[1]*cossin
    
    tuCenter = -1*Au/(2*Auu)
    tvCenter = -1*Av/(2*Avv)
    wCenter = Ao - Auu*tuCenter*tuCenter - Avv*tvCenter*tvCenter
    uCenter = tuCenter*cost - tvCenter*sint
    vCenter = tuCenter*sint + tvCenter*cost
    
    Ru = -1*wCenter/Auu
    Rv = -1*wCenter/Avv
    
    Ru = np.sqrt(np.abs(Ru))*np.sign(Ru)  
    Rv = np.sqrt(np.abs(Rv))*np.sign(Rv)
    
    return np.array([uCenter, vCenter, Ru, Rv, thtarad])

def deleteElement(a1, a2):
    # delete the element of a2 in a1
    a = list(a1)
    for i in a2:
        try:
            a.remove(i)
        except:
            continue
    return np.array(a)

def multiGuass(x, *params):
    # define the mulpiple guass like function
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -(x - ctr)**2/(2*wid**2))
    return y

def multiGuass2(x, *params):
    # define the mulpiple guass like function
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        area = params[i+1]
        wid = params[i+2]
        y = y + area/(wid*np.sqrt(2)*np.pi) * np.exp( -(x - ctr)**2/(2*wid**2))
    return y

def find_repeat(mangaid=None, recommend=True):
    # read the manga repeat_observation.txt to check and return all the repeat
    # observations
    # return the plateifu list, also return the recommend one if recommend is True
    txt = np.loadtxt('data/repeat_observations.txt', skiprows=5, dtype='str')
    target = txt[txt[:, 0] == mangaid]
    repeat_list = []
    if len(target) > 0:
        for i in range(1, 9, 2):
            plate = target[0][i]
            ifu = target[0][i+1]
            if plate != '-1':
                repeat_list.append(plate + '-' + ifu)
        if recommend:
            # pick all the ifu out and find the maximum
            ifu_list = list(map(lambda s: re.match(r"\d+-(\d+)\d", s).group(1), 
                                repeat_list))
            recommend_plateifu = np.array(repeat_list)[np.array(ifu_list)
                                                        == max(ifu_list)]
            # when serveral alternative, the first one is returned
            return repeat_list, recommend_plateifu[0] 
    return repeat_list, None

def dust_corr(Ha, Hb, line_wave, line_flux, E_BV=None, Ha2Hb=3.15, Rv=3.1, quiet=True):
    """ this function using (Ha / Hb) line ration to get the E(B-V) of the 
        reddening curves
    
    Params: 
        Ha, Hb: the flux of the Ha and Hb
        line_wave: the wavelength of correcting line, in um
        line_flux: the flux of the correcting line
        Ha2Hb: the theoretical ratio of Ha and Hb, default to 3.08 for AGN 
               and 2.86 for starforming galaxy
        Rv: effective total obscuration at V, default to 3.1 for Galactic 
            diffused ISM
    return:
        flux: the corrected flux the given emission line
    """
    if E_BV is None:
        Hb = np.ma.masked_less_equal(Hb, 1e-8).filled(np.inf)
        ratio_obs = np.ma.masked_less_equal(Ha/Hb, 1e-8).filled(Ha2Hb)
        # filled the masked value with 9999, to get rid of the warning  message
        E_BV = 1.97*np.log10(ratio_obs / Ha2Hb)
        if not quiet:
            if (E_BV < 0).any():
                print("Unphysical redcorr encountered!")
    E_BV = np.asarray(E_BV)
    E_BV[E_BV < 0] = 0
    # reddening curves using C00 (Calzetti 1999)
    if line_wave > 2.2 or line_wave < 0.12:
        raise ValueError("Wavelenth: {} is out of bound (0.12, 2.2)".format(
                         line_wave))
    elif line_wave < 0.63:
        k_lambda = 2.659*(-2.156+1.509/line_wave-0.198/line_wave**2 \
                   + 0.011/line_wave**3 ) + Rv
    else:
        k_lambda = 2.659*(-1.857+1.040/line_wave)+Rv
    flux_corr = line_flux*10**(0.4*k_lambda*E_BV)
    return np.ma.array(flux_corr, mask=(flux_corr > 20*line_flux))

def curve_parallel(x, d, func, dfunc):
    """generate the parallel data of the given curve

    Params:
        x: the variable
        d: the distance to give curve
        func: the function of the curve
        dfunc: derivative function of the curve
    Return:
        the data parallel to the curve
    """
    xn = x + d*dfunc(x)/np.sqrt(1+dfunc(x)**2)
    yn = func(x) - d/np.sqrt(1+dfunc(x)**2)
    return xn, yn

def curve_dist(x, y, xlim=(-2,-0.1), func=None, dfunc=None, fitfunc=None, p0=None, step=0.01, 
               plot=False):
    """find the distance to the given curve

    Params:
        x, y: the given x coordinates and y coordinates data, 1D array
        func: the curve function
        dfunc: the derivative of func
        p0: the initial value of the fitting parameters
        step: the step interval to generate the paralell curve lines
    Return: distance of x to the curve
    """
    assert len(x) == len(y), "input array should have the length"
    xmin, xmax = x.min(), x.max()
    line_x = np.linspace(*xlim, 100)
    ymin, ymax = y.min(), y.max()
    bins = int(np.rint(max((ymax-ymin)/step, (xmax-xmin)/step)))
    y_fit = np.zeros((bins, len(y)))
    d_interval = np.linspace(step, step*bins, bins)
    
    for i in range(bins): # for different distance
        xn, yn = curve_parallel(line_x, -step*i, func, dfunc)
        popt, pcov = curve_fit(fitfunc, xn, yn, p0=p0)
        y_fit[i,:] = fitfunc(x, *popt)
        # update the initial value for next fitting
        p0 = popt
    distance = abs(y_fit - y)
    min_distance = np.min(distance, axis=0)
    min_d = np.zeros_like(y)
    for i in range(len(y)):
        min_d[i] = d_interval[distance[:,i] <= min_distance[i]][0]
    #plot the func
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(line_x, func(line_x), 'r--')
        im = ax.scatter(x, y ,c=min_d, cmap='RdBu')
        fig.colorbar(im)
    return min_d

def fixmap(fmap, sigma=20, iters=5, filled=1e-8):
    # mask out the unreal value of a given map
    # this function will change the origin array, use with cautious
    if fmap.ndim == 2:
        return sigma_clip(fmap, sigma=sigma, iters=iters).filled(filled)
    elif fmap.ndim == 3:
        if isinstance(1. * filled, float):
            filled = np.repeat(filled, len(fmap))
        assert len(fmap) == len(filled), "filled value and matrix mis-match"
        for k in range(len(fmap)):
            fmap[k] = sigma_clip(fmap[k], sigma=sigma).filled(filled[k])
        return fmap

def bundle_edge(ifudsgn):
    pos = yanny.yanny(filename=package_path+'/data/manga_simbmap_127.par')
    n = (int(ifudsgn)//100 - 1) // 6
    for i in range(1, 8):
        n = n - i
        if n == 0:
            break
    i_low = ((i - 1) * i) * 3 + 1
    i_up = ((i + 1) * i) * 3 + 1
    # print(i_low, i_up)
    pos_ra = pos['SIMBMAP']['raoff'][i_low: i_up]
    pos_dec = pos['SIMBMAP']['decoff'][i_low: i_up]
    pos_ra.append(pos_ra[0])
    pos_dec.append(pos_dec[0])
    return np.array([pos_ra, pos_dec])

def bundle_fibers(ifudsgn):
    pos = yanny.yanny(filename=package_path+'/data/manga_simbmap_127.par')
    n = (int(ifudsgn)//100 - 1) // 6
    for i in range(1, 8):
        n = n - i
        if n == 0:
            break
    i_low = ((i - 1) * i) * 3 + 1
    i_up = ((i + 1) * i) * 3 + 1
    # print(i_low, i_up)
    pos_ra = pos['SIMBMAP']['raoff'][0: i_up]
    pos_dec = pos['SIMBMAP']['decoff'][0: i_up]
    return np.array([pos_ra, pos_dec])

