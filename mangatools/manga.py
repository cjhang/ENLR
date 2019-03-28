#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import matplotlib.backends.backend_pdf as mpdf
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM

from .config import PRODOCTS, DAP_VERSION, DRP_VERSION, SAS
import mangatools as mangatools_package 

DAP_DIR = SAS + '/mangawork/manga/spectro/analysis/' + DRP_VERSION \
              + '/' + DAP_VERSION
MAPS_DIR = DAP_DIR + '/{}-GAU-MILESHC/'.format(PRODOCTS)
DRP_DIR = SAS + '/mangawork/manga/spectro/redux/' + DRP_VERSION
IMAGE_DIR = SAS + "/mangawork/manga/spectro/redux/" + DRP_VERSION \
                + "/allimages/"
DRP = Table.read(DRP_DIR + '/drpall-' + DRP_VERSION + '.fits')
DAP = Table.read(DAP_DIR + '/dapall-' + DRP_VERSION + '-' \
                         + DAP_VERSION + '.fits')

package_path = os.path.dirname(os.path.realpath(mangatools_package.__file__))

class MaNGA(object):
    """Class for a MaNGA observated galaxy"""

    def __init__(self, *arg, mangaid=None):
        """use plateifu and mangaid as the keyword to define a galaxy
        
            *arg: Only the first one is used, as the plateifu
            mangaid: the mangaid of the galaxy

            return MaNGA
        """

        ## Define the global data directories
        self.COSMO = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        self.ESP = 1e-8
        self.C = 299792 # speed of light (km/s)
        ## Pre-read the DRP data

        if arg:
            self.plateifu = arg[0]
            self.drp = DRP[DRP['plateifu']==self.plateifu]
            self.mangaid = self.drp['mangaid'][0]
        elif mangaid:
            self.mangaid = mangaid
            self.drp = DRP[DRP['mangaid']==self.mangaid]
            self.plateifu = self.drp['plateifu'][0]
        else:
            print("No valid plateifu or mangaid detected!")
            raise ValueError
        self.plate, self.ifudsgn = re.split('-', self.plateifu)
        self.mapsfile = MAPS_DIR+self.plate + '/'+self.ifudsgn+'/' \
                        + 'manga-' + self.plate + '-' + self.ifudsgn \
                        + '-MAPS-{}-GAU-MILESHC.fits.gz'.format(PRODOCTS)
        self.datacubefile = MAPS_DIR + '/'+self.plate+'/' \
                            + self.ifudsgn+'/'+'manga-' + self.plate \
                            + '-' + self.ifudsgn \
                            + '-LOGCUBE-{}-GAU-MILESHC.fits.gz'.format(PRODOCTS)
        self.IMAGE_DIR = IMAGE_DIR
        self.ra = self.drp['objra'][0]
        self.dec = self.drp['objdec'][0]
        self.ifura = self.drp['ifura'][0]
        self.ifudec = self.drp['ifudec'][0]
        self.elpetro_r = self.drp['nsa_elpetro_th50_r'][0] 
        self.elpetro_ba = self.drp['nsa_elpetro_ba'][0] 
        self.elpetro_phi = self.drp['nsa_elpetro_phi'][0]
        self.sersic_r = self.drp['nsa_sersic_th50'][0]
        self.sersic_ba = self.drp['nsa_sersic_ba'][0]
        self.sersic_phi = self.drp['nsa_sersic_phi'][0]
        self.sersic_n = self.drp['nsa_sersic_n'][0]
        self.psf = self.drp['rfwhm'][0]
        self.z = self.drp['nsa_z'][0]
        self.c = const.c.to(u.km/u.s).value
        self.d = self.COSMO.comoving_transverse_distance(self.z)
        self.arcsec2kpc = self.COSMO.kpc_comoving_per_arcmin(self.z).to(
                u.kpc/u.arcsec)
        self.repeat, self.alter = self.find_repeat(self)
    
    @property
    def target_range(self):
        # check the target type of primary and secondary
        # currently we still use the 'full fill the field of view' and 'other'
        # by give a resonable value, use with cautious
        mngtarg1 = self.drp['mngtarg1']
        if np.bitwise_and(mngtarg1, 2**10):
            Rrange = 1.5 # 1.5Re Primary
        elif np.bitwise_and(mngtarg1, 2**11):
            Rrange = 2.5 # 2.5Re Secondary
        elif np.bitwise_and(mngtarg1, 2**13):
            Rrange = 2.5 # full fill the field of view
        else:
            Rrange = 1.5 # Other
        return Rrange
    
    def find_repeat(self, recommend=True):
        # read the manga data/repeat_observation.txt to check and return all the 
        # repeat observations
        # return the plateifu list, also return the recommend one if recommend 
        # is True
        mangaid = self.mangaid
        txt = np.loadtxt(package_path + '/data/repeat_observations.txt', skiprows=5, 
                         dtype='str')
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
                ifu_list = list(map(
                        lambda s: re.match(r"\d+-(\d+)\d", s).group(1), 
                        repeat_list))
                recommend_plateifu = np.array(repeat_list)[np.array(ifu_list)
                                                            == max(ifu_list)]
                # when serveral alternative, the first one is returned
                return repeat_list, recommend_plateifu[0] 
        return repeat_list, None

    def image(self, ax=None, showImage=True):
        """return or show rgb-image
        """
        try:
            imagedata = mpimg.imread(IMAGE_DIR + str(self.plate) + '-' + str(self.ifudsgn) + '.png')
        except:
            print("{} image file doesn't exist!".format(self.plateifu))
            imagedata = np.zeros((2,2))

        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.imshow(imagedata)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        # set the 'show=False' when use it as a figure modifier
        if showImage:
            plt.show()

    def showdetail(self, showDRPInfo=True, showDAPInfo=True):
        """ show dap and drp detail of the target
        drp3qual see: https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-5/metadata
        """
        if showDRPInfo:
            drp_qual_dict = dict({0:'GOOD', 1:'RETIRED', 2:'BADDEPTH', 
                4:'SKYSUBBAD', 8:'HIGHSCAT', 16:'BADASTROM', 32:'VARIABLELSF', 
                64:'BADOMEGA', 256:'BADFLUX', 512:'BADPSF', 2**30:'CRITICAL'})
            drp3qual = self.drp['drp3qual'][0]
            print("****** DRP Imformation ******** ")
            print('The drp3qual is: {0}'.format(drp_qual_dict[drp3qual]))
            print('The image information:')
            print('elpetro_r = {0}, b/a = {1}, phi = {2}'.format(
                  self.elpetro_r, self.elpetro_ba, self.elpetro_phi))
            print('The redshift is: {0}'.format(self.z))
            print('The distance is: {0}'.format(self.d))
            print('Re is: {0}'.format(self.arcsec2kpc*self.elpetro_r*u.arcsec))
            # print('LOG10(LOIII) is: {0}'.format(logO3))
        if showDAPInfo:
            dap_qual_dict = dict({0:'GOOD', 1:'FORESTAR', 2:'BADZ', 
                                  4:'LINELESS', 8:'PPXFFAIL',16:'SINGLEBIN', 
                                  2**28:'DRPCRIT', 2**29:'DAPCRIT', 
                                  2**30:'CRITICAL'})
            print("****** DAP Imformation ******** ")
            print('The DAP quality is: {0}'.format(
                    dap_qual_dict[self.maps[0].header['DAPQUAL']]))
    
