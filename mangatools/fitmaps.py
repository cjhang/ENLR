#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The class for reading full spectrum fitting results from datacube.flux_map
"""

import numpy as np
from astropy.io import fits

from astropy.convolution import convolve, Gaussian2DKernel, interpolate_replace_nans
from astropy.stats import sigma_clip

from .maps import Maps

class FitMaps(Maps):
    """class used for reading datacube fitting map
    """
    def __init__(self, plateifu, fitmaps_dir=None, fitmaps_binned_dir=None):
        super().__init__(plateifu)
        self.fitmaps = fits.open(fitmaps_dir + plateifu + '.fits')
        self.fitmaps2 = fits.open(fitmaps_binned_dir + plateifu + '.fits')
        self.agn, *others = self.bptregion()
        
    def line(self, name):
        # acess emission lines
        data = self.fitmaps[name].data
        return data
    
    def line2(self, name):
        # acess emission lines
        data = self.fitmaps2[name].data
        return data
    
    def O3map(self, redcorr=True, smooth=False, fix_outlier=True):
        Ha, Ha_err = self.fitmaps2['Halpha'].data[0], self.fitmaps2['Halpha'].data[1]
        Hb, Hb_err = self.fitmaps2['Hbeta'].data[0], self.fitmaps2['Hbeta'].data[1]
    
        ratio_theory = 3.1
        ratio_obs = Ha / self.fix_zeros(Hb)
        snr_cut = ((Ha / Ha_err) < 3) | ((Hb / Hb_err) < 3)
        ratio_obs[snr_cut] = ratio_theory
        #print("bad red correction pixels:", np.sum(ratio_obs < 3.15))
        #ratio_obs = np.ma.masked_less_equal(ratio_obs, 0).filled(self.Ha2Hb)
        E_BV = 1.97 * np.log10(ratio_obs / ratio_theory)
        E_BV[E_BV < 0] = 0
        E_BV_err = np.sqrt((0.855*Ha_err/self.fix_zeros(Ha))**2 + (0.855*Hb_err/self.fix_zeros(Hb))**2)
        
        if smooth: # smooth the E(B-V) map, due to spaxel were not independant and romove bad pixel
            kernel = Gaussian2DKernel(x_stddev= 0.5 * self.psf/2.355)
            E_BV = convolve(E_BV, kernel, mask=~self.agn, boundary='extend')
            E_BV_err = np.sqrt((0.855*Ha_err/self.fix_zeros(Ha))**2 + (0.855*Hb_err/self.fix_zeros(Hb))**2)
            E_BV_err = convolve(E_BV_err, kernel, mask=~self.agn, boundary='extend')
        
        #if fix_outlier:
        #    corrector = self.fix_outlier(E_BV)
        #    corrector_err = self.fix_outlier(E_BV_err)
            
        k_lambda = 3.52
        corrector = 10**(0.4 * k_lambda * E_BV)
        corrector_err = 0.4 * 10**(0.4 * k_lambda * E_BV) * np.log(10) * k_lambda * E_BV_err
        
        if fix_outlier:
            corrector = self.fix_outlier(corrector)
            corrector_err = self.fix_outlier(corrector_err)
            
        O3 = self.fitmaps['[OIII]5008'].data[0] + self.fitmaps['_[OIII]5008'].data[0]
        O3_err = np.sqrt(self.fitmaps['[OIII]5008'].data[1]**2 + self.fitmaps['_[OIII]5008'].data[1]**2)
        O3_corr = O3 * corrector
        O3_corr_err = np.sqrt((corrector_err * O3)**2 + (O3_err * corrector)**2)
        self.E_BV = E_BV
        self.E_BV_err = E_BV_err
        self.ratio_obs = ratio_obs 
        self.corrector = corrector
        self.corrector_err = corrector_err
        return O3_corr, O3_corr_err


    def fix_zeros(self, arr, filled=True):
        if filled:
            return np.ma.masked_less_equal(arr, 0).filled(np.inf)
        else:
            return np.ma.masked_less_equal(arr, 0)
        
    def fix_outlier(self, arr, interpolate=True, sigma=5, iters=2):
        kernel = Gaussian2DKernel(x_stddev= self.psf/2.355)
        new_array = sigma_clip(arr, sigma=sigma, iters=iters) # remove 5 sigma outlier
        if interpolate:
            new_array = interpolate_replace_nans(new_array.filled(np.nan), kernel)
        return new_array
