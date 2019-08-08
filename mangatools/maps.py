import re
import io
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import matplotlib.backends.backend_pdf as mpdf
from matplotlib.colors import ListedColormap, BoundaryNorm
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.table import Table
from astropy import wcs, coordinates
from astropy.stats import sigma_clip
import requests
import PIL

from .manga import MaNGA
from . import utils

class Maps(MaNGA):
    """ used for emssion line research, based on MaNGA
    """

    def __init__(self, *arg, mangaid=None, mode='G'):
        """
        Args:
            mode: 'G' stands for Gauss fitting data
                  'S' means just make sum of the emission line's flux
        """
        super().__init__(*arg, mangaid=mangaid)
        emline_flux = 'EMLINE_{mode}FLUX'
        emline_ew = 'EMLINE_{mode}EW'
        self.emline_name = emline_flux.format(mode=mode)
        self.emline_ew = emline_ew.format(mode=mode)
        self.maps = fits.open(self.mapsfile)
        self.emline_qual = self.maps[self.emline_name].header['QUALDATA']
        self.emline_err = self.maps[self.emline_name].header['ERRDATA']
        self.emline_ew_qual = self.maps[self.emline_ew].header['QUALDATA']
        self.emline_ew_err = self.maps[self.emline_ew].header['ERRDATA']
        self.emline_dict = utils.emline_dictionary(self.maps, self.emline_name)
        self.naxis1 = self.maps[1].header['NAXIS1'] 
        self.naxis2 = self.maps[1].header['NAXIS2']
        self.mask = self.maps['BINID'].data[0] < 0 # continuum binid
        self.ra_map, self.dec_map = np.ma.array(self.maps['SPX_SKYCOO'].data,
                                                mask=[self.mask, self.mask])
        #self.x, self.y= np.flip(x), y
        # Mean g-band-weighted signal-to-noise bin
        self.binmatrix = self.maps['BIN_SNR'].data
        self.header = self.maps['SPX_MFLUX'].header
        self.wcs_map = wcs.WCS(self.maps['SPX_MFLUX'].header)

    def skycoord_to_pixel(self, coord):
        xy = wcs.utils.skycoord_to_pixel(coord, self.wcs_map)
        return np.round(xy).astype(int)

    def redcorrmap(self, Ha2Hb=3.1, smooth=True, plot=True, showImage=True, 
                   mini=False, ax=None, showColorbar=True, fs=10, **kwargs):
        # return the E_BV map of the galaxy
        Ha = self.line('Ha-6564', snr=3, redcorr=False)
        Hb = self.line('Hb-4862', snr=3, redcorr=False)
        ratio_obs = np.ma.array(Ha/Hb)
        # filled the masked value with 9999, to get rid of the warning  message
        ratio_obs2 = np.ma.array(ratio_obs.filled(9999), mask=ratio_obs.mask)
        E_BV = 1.97*np.log10(ratio_obs2 / Ha2Hb)
        E_BV[E_BV < 0] = 0 # remove the unphysical correction
        #E_BV.filled(0)
        if smooth: # smooth the E(B-V) map, due to spaxel were not idependant
            kernel = Gaussian2DKernel(x_stddev=0.5*self.psf/2.355)
            E_BV = np.ma.array(convolve(E_BV, kernel), mask=ratio_obs.mask)
        if plot:
            if ax == None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            x, y = self.maps['SPX_SKYCOO'].data
            x = np.flip(x)
            im = ax.pcolormesh(x, y, E_BV, **kwargs)
            if mini:
                ax.tick_params(axis='both', which='major', labelsize=fs)
                # plot.tick_params(axis='both', which='minor', labelsize=8)
                # ax.set_xticklabels([])
                # ax.set_yticklabels([])
            else:
                ax.set_title('E(B-V)', fontsize=fs)
                ax.set_ylabel('arcsec', fontsize=fs)
                ax.set_xlabel('arcsec', fontsize=fs)
            if showColorbar:
                if mini:
                    cbar = plt.colorbar(im, ax=ax, pad=0.0)
                else: 
                    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=fs)
            if showImage:
                return fig
        else:
            return E_BV
   
    def image_quary(self, ax=None, showImage=True, mini=True, scale=0.1,
                    width=640, height=640, band='R', opt='G', lw = 1, fs=8, 
                    mag_range=(0, 21), show_bundle=True, showid=True, 
                    show_fibers=False, highlight=None):
        base_url = 'http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra}&dec={dec}&scale={scale}&width={width}&height={height}&opt={opt}&query={band} {mag_range}'
        url = base_url.format(ra=self.ifura, dec=self.ifudec, scale=scale,
                              width=width, height=height, band=band, opt=opt,
                              mag_range=mag_range)
        r = requests.get(url)
        # print(r.url)
        if r.status_code != 200:
            raise ValueError("request error! url:{}".format(url))
        im = PIL.Image.open(io.BytesIO(r.content))
        ax1, ax2 = im.size
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.imshow(im)
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if show_bundle:
            bundle_edge = utils.bundle_edge(self.ifudsgn)
            ax.plot(ax1 / 2 + 1 / scale * bundle_edge[0], 
                    ax2 / 2 + 1 / scale * bundle_edge[1], 
                    color='magenta', lw=lw, alpha=1)
        if show_fibers:
            fibers_pos = utils.bundle_fibers(self.ifudsgn)
            ax.plot(ax1 / 2 + 1 / scale * fibers_pos[0], 
                    ax2 / 2 + 1 / scale * fibers_pos[1], 
                    'o', color='white', lw=lw, ms=1/scale-6, 
                    fillstyle='none', alpha=0.2)
            if highlight is not None:
                ax.plot(ax1 / 2 + 1 / scale * fibers_pos[0][highlight], 
                        ax2 / 2 + 1 / scale * fibers_pos[1][highlight], 
                        'o', color='red', lw=lw+2, ms=1/scale-6, 
                        fillstyle='none', alpha=0.5)

        if showid:
            ax.text(width*0.5, height*0.1, 'plateifu: '+self.plateifu, color='white', fontsize=fs, ha='center')
        if showImage:
            if ax == None:
                return fig
            #plt.show()
    
    def image_quary2(self, ax=None, showImage=True, mini=True, scale=0.1,
                    width=640, height=640, band='R', opt='G',
                    mag_range=(0, 21), show_bundle=True):
        base_url = 'http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra}&dec={dec}&scale={scale}&width={width}&height={height}&opt={opt}&query={band} {mag_range}'
        url = base_url.format(ra=self.ifura, dec=self.ifudec, scale=scale,
                              width=width, height=height, band=band, opt=opt,
                              mag_range=mag_range)
        r = requests.get(url)
        # print(r.url)
        if r.status_code != 200:
            raise ValueError("request error! url:{}".format(url))
        im = PIL.Image.open(io.BytesIO(r.content))
        ax1, ax2 = im.size
        self.ax1, self.ax2 = ax1, ax2
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.imshow(im)
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if show_bundle:
            bundle_edge = utils.bundle_edge(self.ifudsgn)
            ax.plot(ax1 / 2 + 1 / scale * bundle_edge[0], 
                    ax2 / 2 + 1 / scale * bundle_edge[1], 
                    color='magenta', lw=1)
        if showImage:
            if ax == None:
                return fig
            #plt.show()

    def showpoint_image(self, xy, xy0, ax=None, showImage=True, radius=3, 
                        mini=False, scale=0.5):
        xy = np.atleast_1d(xy) # the xy is based on scale = 0.5
        xy0 = np.atleast_1d(xy0)
        xy_center = np.array([self.header['CRPIX1'], self.header['CRPIX2']])
        #xy_new = (xy - xy_center) * 0.5 / scale + xy0
        # print(xy)
        # print(xy0)
        # if len(x,y)
        x_new = (xy[:, 0] - xy_center[0]) * 0.5 / scale + xy0[0]
        y_new = - (xy[:, 1] - xy_center[1]) * 0.5 / scale + xy0[1]
        xy_new = np.array(list(zip(x_new, y_new)))
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        from matplotlib.patches import Circle
        if xy.ndim == 1:
            circle = Circle(xy_new, radius=radius/scale, zorder=2, lw=1,
                            label=xy_new,
                            edgecolor='magenta', fill=False)
            ax.add_artist(circle)
        elif xy.ndim > 1:
            for i in range(len(xy)):
                circle = Circle(xy_new[i], radius=radius/scale, zorder=2, 
                                label=xy_new,
                                edgecolor='C{}'.format(i), 
                                fill=False)
                ax.add_artist(circle)
        else:
            raise ValueError("Invaid coordinates!")
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if showImage:
            if ax == None:
                return fig
    
    def showpoint_image2(self, ax=None, showImage=True, radius=3, 
                        mini=False, scale=0.5, show_Re=True, show_id=True,
                        Re=None, ba=None, phi=None, color='k'):
        xy_center = np.array([self.ax1//2, self.ax2//2])
        #xy_new = (xy - xy_center) * 0.5 / scale + xy0
        # print(xy)
        # print(xy0)
        # if len(x,y)
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if show_id:
            ax.set_title(self.plateifu)
        if show_Re:
            from matplotlib.patches import Ellipse
            ellipse = Ellipse(xy_center, Re * ba, Re, angle=phi, zorder=2, 
                              edgecolor=color, fill=False)
            ax.add_artist(ellipse)
        if showImage:
            if ax == None:
                return fig

    def showpoint(self, xy, linename='OIII-5008', flux='total', ax=None, 
                  show_aperture=True, showImage=True, radius=5, mini=False,
                  show_Re=False, Re=None, ba=None, phi=None, scale=1):
        xy = np.atleast_1d(xy)
        x_coor, y_coor = self.maps['SPX_SKYCOO'].data
        x_coor = np.flip(x_coor)

        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if flux == 'total':
            #ax.pcolormesh(x_coor, y_coor, flux)
            flux_tot = self.maps['SPX_MFLUX'].data #g-band-weighted mean flux
            ax.imshow(flux_tot, origin='lower')
        elif flux:
            ax.imshow(flux, origin='lower')
        if show_aperture:
            from matplotlib.patches import Circle
            if xy.ndim == 1:
                circle = Circle(xy, radius=radius, zorder=2, lw=1,
                                edgecolor='magenta', fill=False)
                ax.add_artist(circle)
            elif xy.ndim > 1:
                for i in range(len(xy)):
                    circle = Circle(xy[i], radius=radius, zorder=2, 
                                    edgecolor='C{}'.format(i), 
                                    fill=False)
                    ax.add_artist(circle)
            else:
                raise ValueError("Invaid coordinates!")
        
        if show_Re:
            from matplotlib.patches import Ellipse
            ellipse = Ellipse(xy, Re * ba, Re, angle=phi, zorder=2, 
                              edgecolor='k', fill=False)
            ax.add_artist(ellipse)
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if showImage:
            if ax == None:
                return fig

    def image(self, ax=None, showImage=True, mini=False, showid=True, fs=8):
        """return or show rgb-image
        """
        try:
            imagedata = mpimg.imread(self.image_file)
        except:
            print("{} image file doesn't exist!".format(self.plateifu))
            imagedata = np.zeros((2,2))
        x, y = self.maps['SPX_SKYCOO'].data
        extent = [x.min(), x.max(), y.min(), y.max()]
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.imshow(imagedata, origin='upper', extent=extent)
        if showid:
            ax.text(x.mean(), y.max()*0.8, 'MaNGA ID: '+self.mangaid, color='white', fontsize=fs, ha='center')
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        # set the 'show=False' when use it as a figure modifier
        if showImage:
            return fig
            #plt.show()

    def line(self, name, mask=True, snr=3, redcorr=False, average=True, 
             select_region=None, return_err=False, plot=False, ax=None, 
             showImage=True, show_psf=False, mini=False, fs=10, 
             showColorbar=True, show_bundle=False, **kwargs):
        '''return the masked emission line

        Params:
            average: using a give E_BV to do red correction
        '''
        # get the data
        emline_flux = self.maps[
                self.emline_name].data[self.emline_dict[name]]
        emline_flux_ivar = self.maps[
                self.emline_err].data[self.emline_dict[name]]
        emline_flux_err = np.sqrt(1/(emline_flux_ivar+1e-8))
        if snr:
            emline_flux_snr = emline_flux*np.sqrt(emline_flux_ivar)
            snr_mask = emline_flux_snr < snr
        else:
            snr_mask = None
        if mask:
            emline_flux_mask = self.maps[self.emline_qual].data[
                    self.emline_dict[name]] > 0
        else:
            emline_flux_mask = None
        final_mask = np.logical_or(snr_mask, emline_flux_mask)
        if redcorr:
            #Ha = self.maps[self.emline_name].data[self.emline_dict['Ha-6564']]
            #Hb = self.maps[self.emline_name].data[self.emline_dict['Hb-4862']]
            Ha = self.line('Ha-6564', snr=3, redcorr=False)
            Hb = self.line('Hb-4862', snr=3, redcorr=False)
            line_wave = re.match(r"\w+-(\d+)", name)[1]
            if average:
                Ha2Hb=3.15
                ratio_obs = np.ma.array(Ha/Hb)
                # filled the masked value, to get rid of the warning  message
                ratio_obs2 = np.ma.array(ratio_obs.filled(9999), 
                                         mask=ratio_obs.mask)
                E_BV = 1.97*np.log10(ratio_obs2 / Ha2Hb)
                E_BV = np.mean(E_BV[select_region])
                emline_flux = utils.dust_corr(0, 0, float(line_wave)/10000, 
                                              emline_flux, E_BV=E_BV)
            else:
                emline_flux = utils.dust_corr(Ha, Hb, float(line_wave)/10000, 
                                              emline_flux)
        if plot:
            if ax == None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            flux = self.line(name, **kwargs)
            x, y = self.maps['SPX_SKYCOO'].data
            x = np.flip(x)
            # x_coor, y_coor = self.maps['SPX_SKYCOO'].data
            extent = [x.min(), x.max(), y.min(), y.max()]
            # im = ax.imshow(flux, extent=extent, **kwargs)
            im = ax.pcolormesh(x, y, flux, cmap='viridis', **kwargs)
            if mini:
                ax.tick_params(axis='both', which='major', labelsize=fs)
                # plot.tick_params(axis='both', which='minor', labelsize=8)
                # ax.set_xticklabels([])
                # ax.set_yticklabels([])
            else:
                ax.set_title('{} Flux'.format(name), fontsize=fs)
                ax.set_ylabel('arcsec', fontsize=fs)
                ax.set_xlabel('arcsec', fontsize=fs)
            cs = ax.contour(x, y, flux, levels=[1, 3,10,100], colors='r', 
                       extent=extent)
            plt.clabel(cs, inline=1, fontsize=5, fmt='%1.1f')
            if showColorbar:
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=fs)
            if show_psf:
                from matplotlib.patches import Circle
                ellipse = Circle((x.mean(), y.mean()), radius=self.psf, 
                                 zorder=2, edgecolor='lightpink', fill=False)
                ax.add_artist(ellipse)
            if show_bundle:
                self.bundle_edge(ax)
            if showImage:
                return fig
                #plt.show(fig)
        else:
            if return_err: # return the error^2
                return np.ma.MaskedArray(emline_flux_err, mask=final_mask, 
                                         fill_value=99999)
            return np.ma.MaskedArray(emline_flux, mask=final_mask, 
                                     fill_value=99999)
    
    def lines(self, **kwargs):
        # get the dictionary contains all the lines
        x, y = self.maps['SPX_SKYCOO'].data
        lines_dict = {}
        for emline in self.emline_dict.keys():
            lines_dict[emline] = Maps.line(self, emline, **kwargs)
        return lines_dict

    def line_ew(self, name, mask=True, snr=None, plot=False, ax=None, 
                showImage=True, mini=False, fs=10, showColorbar=True, 
		show_bundle=False, **kwargs):
        emline_ew = self.maps[self.emline_ew].data[self.emline_dict[name]]
        emline_ew_ivar = self.maps[self.emline_ew_err].data[self.emline_dict[name]]
        if snr:
            emline_ew_snr = emline_ew*np.sqrt(emline_ew_ivar)
            snr_mask = emline_ew_snr < snr
        else:
            snr_mask = None
        if mask:
            emline_ew_mask = self.maps[self.emline_ew_qual].data[self.emline_dict[name]] > 0
        else:
            emline_ew_mask = None
        final_mask = np.logical_or(snr_mask, emline_ew_mask)
        if plot:
            if ax == None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            ew = self.line_ew(name, snr=snr)
            x, y = self.maps['SPX_SKYCOO'].data
            x = np.flip(x)
            extent = [x.min(), x.max(), y.min(), y.max()]
            # im = ax.imshow(ew, extent=extent, **kwargs)
            im = ax.pcolormesh(x, y, ew, **kwargs)
            cs = ax.contour(x, y, ew, levels=[3,10,100], colors='r', 
                       extent=extent, linewidths=0.5)
            plt.clabel(cs, inline=1, fontsize=5, fmt='%1.1f')
            if mini:
                ax.set_xticklabels([])
                ax.set_yticklabels([])
            else:
                ax.set_xlabel('arcsec')
                ax.set_ylabel('arcsec')
                ax.set_title('{} EW'.format(name))
            if showColorbar:
                cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=fs)
            if show_bundle:
                self.bundle_edge(ax)
            if showImage:
                return fig
                #plt.show()
        else:
            return np.ma.MaskedArray(emline_ew, mask=final_mask)
        
    def vfield(self, component='stellar', name=None, plot=False, ax=None, 
		showImage=True, show_bundle=False, mini=False, fs=10):
        # Create a masked array with the stellar velocity data
        
        if component == 'stellar':
            mask_ext = self.maps['STELLAR_VEL'].header['QUALDATA']
            stellar_vfield = self.maps['STELLAR_VEL'].data
            mask = self.maps[mask_ext].data > 0
            stellar_vfield = np.ma.MaskedArray(stellar_vfield, mask=mask)
            stellar_vfield = sigma_clip(stellar_vfield)
            vfield = stellar_vfield
        if component == 'gas':
            if name is None:
                raise ValueError('Gas name is required!')
            mask_ext = self.maps['EMLINE_GVEL'].header['QUALDATA']
            gas_vfield = self.maps['EMLINE_GVEL'].data[self.emline_dict[name]]
            vfield = gas_vfield
        # Show the stellar velocity field
        if plot:
            if ax == None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            x, y = self.maps['SPX_SKYCOO'].data
            x = np.flip(x)
            im = ax.pcolormesh(x, y, vfield, cmap='RdBu_r')
            if mini:
                ax.set_xticklabels([])
                ax.set_yticklabels([])
            else:
                ax.set_title('velocity')
                ax.set_xlabel('arcsec')
                ax.set_ylabel('arcsec')
            # fig.colorbar(im, ax=ax)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=fs)
            if show_bundle:
                self.bundle_edge(ax)
            if showImage:
                return fig
                #plt.show()
        return vfield

    def sigmafield(self, ax=None, showImage=True, mini=False, fs=10):
        # Create a masked array with the stellar velocity data
        mask_ext = self.maps['STELLAR_SIGMA'].header['QUALDATA']
        stellar_vfield = self.maps['STELLAR_SIGMA'].data
        mask = self.maps[mask_ext].data > 0
        stellar_vfield = np.ma.MaskedArray(stellar_vfield, mask=mask)
        stellar_vfield = sigma_clip(stellar_vfield)
        # Show the stellar velocity field
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        x, y = self.maps['SPX_SKYCOO'].data
        x = np.flip(x)
        im = ax.pcolormesh(x, y, stellar_vfield, cmap='RdBu_r')
        if mini:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_title('velocity')
            ax.set_xlabel('arcsec')
            ax.set_ylabel('arcsec')
        # fig.colorbar(im, ax=ax)
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=fs)
        if showImage:
            return fig
            #plt.show()

    def bptregion(self, snr=3, mask=True, strickMode=False, mode='N2'):
        x_coor, y_coor = self.maps['SPX_SKYCOO'].data
        # self.image(fig=fig, index=index, showImage=False)
        lines = self.lines(mask=mask, snr=snr)
        if mode == 'N2':
            ## plot O3-N2
            x0 = lines['NII-6585']/lines['Ha-6564']
            y0 = lines['OIII-5008']/lines['Hb-4862']
        elif mode == 'S2':
            x0 = (lines['SII-6718']+lines['SII-6732'])/lines['Ha-6564']
            y0 = lines['OIII-5008']/lines['Hb-4862']
        else:
            raise ValueError
        # filled the masked data with 9999, to avoid the warning massage
        x = np.ma.array(np.log10(x0.filled(9999)), mask=x0.mask)
        y = np.ma.array(np.log10(y0.filled(9999)), mask=y0.mask)
        self.bptx, self.bpty = x, y #convinent for external usage
        return utils.bptregion(x, y, mode=mode)

    def bpt(self, ax=None, showImage=True, show_model=True, mini=False, fs=10, 
            showTitle=True, show_schawinski=False, showLineName=False, **kwargs):
        """1D version of BPT, plot all the spaxels on the BPT diagram
        """
        agn, cp, sf, *others = self.bptregion(**kwargs)

        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if show_model:
            # Ke01 line [Kewly et al (2006)]
            x1 = np.linspace(-1.5, 0.2, 100)
            y_ke01 = 0.61/(x1-0.47)+1.19
            # Ka03 line [Kauffmann et al.(2003)]
            x2 = np.linspace(-1.5, -.2, 100)
            y_ka03 = 0.61/(x2-0.05)+1.3
            # Schawinski2007 imperical line
            x3 = np.linspace(-0.18, 0.7, 100)
            y_schawinski = 1.05*x3+0.45
            ax.plot(x1, y_ke01, '--', zorder=3, color='C0')
            ax.plot(x2, y_ka03, '-', zorder=3, color='C1')
            if showLineName:
                ax.text(-0.1, -0.4, 'ke01', rotation=-70, fontsize=fs, color='C0')
                ax.text(-0.7, 0.2, 'ka03', rotation=-55, fontsize=fs, color='C1')
            if show_schawinski:
                ax.plot(x3, y_schawinski, 'r-.', zorder=3)
            #ax.plot(x1, y_ke01, 'r', zorder=1)
            #ax.plot(x2, y_ka03, 'k--', zorder=2)
            #ax.plot(x3, y_schawinski, 'k--', zorder=3)
            ax.set_xlim((-1.4, 0.8))
            ax.set_ylim((-1., 1.4))
            #ax.set_xlim(-1.25, .5)
            #ax.set_ylim(-0.6, 1.3)
        ax.plot(self.bptx[agn].flatten(), self.bpty[agn].flatten(), 'ro', ms=2, lw=1)
        ax.plot(self.bptx[cp].flatten(), self.bpty[cp].flatten(), 'o',
                color='lightgreen', ms=2, lw=2)
        ax.plot(self.bptx[sf].flatten(), self.bpty[sf].flatten(), 'bo', ms=2, lw=1)
        if mini:
            ax.tick_params(axis='both', which='major', labelsize=fs)
            # ax.set_xticklabels([])
            # ax.set_yticklabels([])
        else:
            ax.tick_params(axis='both', which='major', labelsize=fs)
            ax.tick_params(axis='both', which='minor', labelsize=fs-1)
            ax.set_xlabel(r'LOG([NII]/H$\alpha)$', fontsize=fs)
            ax.set_ylabel(r'LOG([OIII]/H$\beta)$', fontsize=fs)
            if showTitle:
                ax.set_title('BPT', fontsize=fs)
        if showImage:
            return fig
            #plt.show(fig)

    def bpt2d(self, snr=3, mask=True, mode='N2', strickMode=False, ax=None, 
              showImage=True, show_Re=False, Re=1, show_psf=False, fs=8, 
              showTitle=True, mini=False, show_aperture=False, R_ape=3,
              show_bundle=False):
        """
        Args:
            show_half_Re, show_psf all refer to the r band
        """
        #x_coor, y_coor = self.maps['SPX_SKYCOO'].data
        xx = (np.arange(self.naxis1) - self.naxis1 * 0.5) * 0.5
        yy = (np.arange(self.naxis2) - self.naxis2 * 0.5) * 0.5
        x_coor, y_coor = np.meshgrid(xx, yy)
        # x_coor = np.flip(x_coor)
        # self.image(fig=fig, index=index, showImage=False)
        region_type = np.full(x_coor.shape, np.nan)
        lines = self.lines(mask=mask, snr=snr)
        if mode == 'N2':
            ## plot O3-N2
            x0 = lines['NII-6585']/lines['Ha-6564']
            y0 = lines['OIII-5008']/lines['Hb-4862']
        elif mode == 'S2':
            x0 = (lines['SII-6718']+lines['SII-6732'])/lines['Ha-6564']
            y0 = lines['OIII-5008']/lines['Hb-4862']
        else:
            raise ValueError
        # filled the masked data with 9999, to avoid the warning massage
        x = np.ma.array(np.log10(x0.filled(9999)), mask=x0.mask)
        y = np.ma.array(np.log10(y0.filled(9999)), mask=y0.mask)
        if mode == 'N2':
            if strickMode:
                region_color = ['red','lightgreen','blue', 'darkred', 'orange']
                region_name = ['AGN', 'Comp', 'HII', 'Seyfert', 'LINER']
                
                AGN, CP, SF, seyfert, liner= utils.bptregion(x, y, mode='N2')
                region_type[AGN] = 1
                region_type[CP] = 2
                region_type[SF] = 3
                region_type[seyfert] = 4
                region_type[liner] = 5
                bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5] # set color for imshow
            else:
                region_color = ['red','lightgreen','blue'] 
                region_name = ['AGN', 'Comp', 'HII']
                AGN, CP, SF, *not_need= utils.bptregion(x, y, mode='N2')
                region_type[AGN] = 1
                region_type[CP] = 2
                region_type[SF] = 3
                bounds = [0.5, 1.5, 2.5, 3.5] # set color for imshow
        elif mode == 'S2':
            region_color = ['red','orange','blue'] 
            region_name = ['Seyfert', 'LIER', 'HII']
            seyfert, lier, sf = utils.bptregion(x, y, mode='S2')
            region_type[seyfert] = 1
            region_type[lier] = 2
            region_type[sf] = 3
            bounds = [0.5, 1.5, 2.5, 3.5] # set color for imshow

        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        from matplotlib import colors
        import matplotlib.patches as mpatches
        cmap = colors.ListedColormap(region_color)
        norm = colors.BoundaryNorm(bounds, cmap.N)
        # the map data was up-down inverse as the optical image from sdss
        ax.pcolormesh(x_coor, y_coor, region_type, cmap=cmap, norm=norm)
        if mini:
            ax.tick_params(axis='both', which='major', labelsize=fs)
            # ax.set_xticklabels([])
            # ax.set_yticklabels([])
        else:
            ax.tick_params(axis='both', which='major', labelsize=fs)
            ax.tick_params(axis='both', which='minor', labelsize=fs-1)
            ax.set_xlabel('arcsec', fontsize=fs)
            ax.set_ylabel('arcsec', fontsize=fs)
            if showTitle:
                ax.set_title('BPT', fontsize=fs)
        # fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        patches = [mpatches.Patch(color=region_color[i], 
                   label=region_name[i] ) for i in range(len(region_name))]
        ax.legend(handles=patches, loc='upper right', borderaxespad=0., 
                  fontsize=fs)
        if show_Re:
            from matplotlib.patches import Ellipse
            ellipse_r = Re*self.elpetro_r
            ellipse = Ellipse((0, 0), 
                               ellipse_r*self.elpetro_ba*2, ellipse_r*2, 
                               angle=self.elpetro_phi, zorder=2, 
                               edgecolor='orange', fill=False)
            ax.add_artist(ellipse)
        if show_psf:
            from matplotlib.patches import Circle
            circle = Circle((0, 0), radius=self.psf/2, 
                             zorder=2, edgecolor='darkviolet', fill=False)
            ax.add_artist(circle)
        if show_aperture:
            from matplotlib.patches import Circle
            circle = Circle((0, 0), radius=R_ape, 
                             zorder=2, edgecolor='navy', fill=False)
            ax.add_artist(circle)
        if show_bundle:
            self.bundle_edge(ax)
        if showImage:
            return fig
            #plt.show(fig)

    def eshape(self, plot=False, ax=None, fs=10, lw=1, ms=2 ,Re=None, 
               useHaEW=False):
        """calculate the extended shape
        
        the morphology of ENLR was distinguished by human before, here try to 
        classified it by fourier series following (Zhicheng He 2018)
        """
        if Re is None:
            Re = self.target_range
        frac_list = []
        phi_list = []
        r, rr, phi = self.maps['SPX_ELLCOO'].data
        agn, cp, sf, seyfert, liner, *other = self.bptregion()
        if useHaEW:
            HaEW_region = self.line_ew('Ha-6564') > 3
            agn = agn & HaEW_region
        for i in range(36):
            phi_min = 10*i
            phi_max = phi_min+10
            phi_list.append(phi_min+5)
            section_bin = utils.sector_binning(r, phi, r_min=0, 
                    r_max=Re*self.elpetro_r, phi_min=phi_min, phi_max=phi_max)
            #return section_bin, l.maps['BINID'].data[0]
            total_pixel = np.sum((self.maps['BINID'].data[0]>0) & section_bin)
            if total_pixel == 0:
                #print("No valid area of {}".format(plateifu))
                #return 0, 0
                frac = 0
            else:
                #frac = np.sum(agn & section_bin)/np.sum(l.maps['BINID'].data[0] & section_bin)
                frac = np.sum(agn & section_bin)/total_pixel
            frac_list.append(frac)
        #a0 = np.sum(frac_list)/18
        a0 = np.mean(frac_list)
        #print(np.array(frac_list).shape, np.array(phi_list).shape)
        a1 = np.sum(np.array(frac_list)*np.cos(1*np.array(phi_list)/180*np.pi))/18
        b1 = np.sum(np.array(frac_list)*np.sin(1*np.array(phi_list)/180*np.pi))/18
        a2 = np.sum(np.array(frac_list)*np.cos(2*np.array(phi_list)/180*np.pi))/18
        b2 = np.sum(np.array(frac_list)*np.sin(2*np.array(phi_list)/180*np.pi))/18
        # plot the fitting data
        if plot:
            if ax == None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            ax.step(phi_list, frac_list, lw=lw)
            ax.tick_params(axis='both', which='major', labelsize=fs)
            ax.tick_params(axis='both', which='minor', labelsize=fs-1)
            ax.set_xlabel('phi', fontsize=fs)
            ax.set_ylabel('fraction', fontsize=fs)
            yfit0 = a0
            # yfit1 = f(phi_list, a1, b1, 1)+a0
            # yfit2 = f(phi_list, a2, b2, 2)+a0
            theta = np.array(phi_list) / 180 * np.pi
            yfit1 = a1 * np.cos(theta) + b1 * np.sin(theta) + a0
            yfit2 = a2 * np.cos(2 * theta) + b2 * np.sin(2 * theta) + a0
            ax.plot(phi_list, np.zeros_like(phi_list)+yfit0, ':',label=r'$a_0$', lw=lw, ms=ms)
            ax.plot(phi_list, yfit1, '--',label=r'$C_{m=1} + a_0$', lw=lw, ms=ms)
            ax.plot(phi_list, yfit2, '-.', label=r'$C_{m=2} + a_0$', lw=lw, ms=ms)
            ax.legend(loc='upper right', fontsize=fs)
            ax.text(10, np.max(frac_list)+0.22, 
                    "\n$\lambda_1$={:.2f} \n$\lambda_2$={:.2f}".format(np.sqrt(a1**2+b1**2)/a0, 
                                                          np.sqrt(a2**2+b2**2)/a0),
                    horizontalalignment='left', verticalalignment='top', fontsize=fs+2)
            ax.set_ylim(np.min(frac_list)-0.1, np.max(frac_list)+0.2)
        # print('a0 is {}'.format(a0))
        return np.sqrt(a1**2+b1**2)/a0, np.sqrt(a2**2+b2**2)/a0

    def plot(self, lines=['OIII-5008',], lines_ew=['Ha-6564',], 
             cols=4, showImage=True, mode='N2'):
        """fast plot
        """
        subplot_nums = len(lines) + len(lines_ew) + 3
        if (subplot_nums <= cols) or (cols == 1):
            fig, ax_array = plt.subplots(1, subplot_nums, 
                                         sharex=True, sharey=True, 
                                         figsize=(3*subplot_nums, 3))
            indx = 0
            #for image
            self.image(ax=ax_array[indx], showImage=False)
            indx += 1
            # for BPT
            self.bpt2d(ax=ax_array[indx], snr=3, strickMode=False, mode=mode,
                       showImage=False, show_aperture=True, show_psf=True)
            indx += 1
            # for velocity field
            self.vfield(ax=ax_array[indx], showImage=False, plot=True)
            indx += 1
            for i in range(len(lines)):
                self.line(lines[i], plot=True, ax=ax_array[indx], 
                          showImage=False, snr=3, redcorr=False)
                indx += 1
            for i in range(len(lines_ew)):
                self.line_ew(lines_ew[i], plot=True, ax=ax_array[indx], 
                          showImage=False)
                indx += 1
        # lines_rows = int(np.ceil(len(lines)/cols))
        # lines_ew_rows = int(np.ceil(len(lines_ew)/cols))
        # subplot_rows = lines_rows + lines_ew_rows + 1
        else:
            def gen(indx, indy): # generate the column and row
                if indx > cols -2:
                    indx = 0
                    indy += 1
                else:
                    indx += 1
                return indx, indy
            subplot_rows = int(np.ceil(subplot_nums/cols))
            fig, ax_array = plt.subplots(subplot_rows, cols, sharex=True, 
                                sharey=True, figsize=(4*cols, 3*subplot_rows))
            indx, indy = 0, 0
            self.image(ax=ax_array[indy, indx], showImage=False)
            indx, indy = gen(indx, indy)
            self.bpt2d(ax=ax_array[indy, indx], snr=3, mode=mode, 
                    strickMode=False, showImage=False, show_Re=0.5, show_psf=True)
            indx, indy = gen(indx, indy)
            self.vfield(ax=ax_array[indy,indx], showImage=False)
            indx, indy = gen(indx, indy)
            for i in range(len(lines)):
                self.line(lines[i], plot=True, ax=ax_array[indy, indx], 
                          showImage=False, snr=3, redcorr=True)
                indx, indy = gen(indx, indy)
            for i in range(len(lines_ew)):
                self.line_ew(lines_ew[i], plot=True, ax=ax_array[indy, indx], 
                             showImage=False)
                indx, indy = gen(indx, indy)
        plt.suptitle(self.plateifu)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.8, bottom=0.2)
        if showImage:
            plt.show(fig)
        else:
            plt.close(fig)
        return fig

    def bundle_edge(self, ax):
        """show the shadow of bundle coverage"""
        x_coor, y_coor = self.maps['SPX_SKYCOO'].data
        x_coor = np.flip(x_coor)
        with fits.open(self.logcube_file) as f:
            bundle_mask = f['MASK'].data[0, :, :]
        cmap = ListedColormap(['#cfcfcf', '#ffffff'])
        bounds = [-0.5, 0.5, 1.5]
        norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
        ax.pcolormesh(x_coor, y_coor, bundle_mask == 1027, cmap=cmap, norm=norm,
                zorder=-1)

