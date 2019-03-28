import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.table import Table
from astropy import wcs, coordinates

from .ppxf import ppxf, ppxf_util, miles_util, capfit
from vorbin.voronoi_2d_binning import voronoi_2d_binning

from .manga import MaNGA
from .maps import Maps
from .spaxel import Spaxel

class DatacubeMask():
    """Used for generate the mask mastrix for datacube
    """
    # manga dap spectrum mask, see https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-6/DAPMetaData#MANGA_DAPSPECMASK
    MANGA_DAPSPECMASK = dict({'IGNORED':1, 'FORESTAR':2, 'FLUXINVALID':4,
                             'IVARINVALID':8, 'ARTIFACT':16, 'FITIGNORED':32,
                             'FITFAILED':64, 'ELIGNORED':128, 'ELFAILED':256})
    def flagged(self, maskId, maskFlags):
        # try to write a function about the mangadap's DAPCubeBitMask
        # not finish yet
        # flag_masks = list(map(lambda name: self.MANGA_DAPSPECMASK[name], maskFlags))
        mask = np.zeros_like(maskId)
        for flag in maskFlags:
            mask |= np.bitwise_and(maskId, self.MANGA_DAPSPECMASK[flag])
            # mask = np.logical_or(mask, 
                    # np.bitwise_and(maskId, self.MANGA_DAPSPECMASK[flag]))
        return mask

class Datacube(MaNGA):
    """
    test galaxies can be 8141-1901(AGN), 8549-12701, 7815-6104, 7495-1902, 8547-12701
                                         8249-3704, 8597-3703
    """

    def __init__(self, plateifu, mask=True, load_fit=True, redcorr=True):
        super().__init__(plateifu)
        self.mask = mask
        self.load_fit = load_fit
        self.bm = DatacubeMask()
        self.maps = Maps(plateifu)
        self.datacube = fits.open(self.datacubefile)
        self.hdu = self.datacube['FLUX']
        self.wcs_cude = wcs.WCS(self.hdu.header)
        self.wave = self.datacube['WAVE'].data
        self.wave_rest = self.wave/(1+self.z)
        self.flux = np.ma.array(
                self.datacube['FLUX'].data, 
                mask=self.bm.flagged(self.datacube['MASK'].data, \
                ['IGNORED', 'FLUXINVALID', 'ARTIFACT']))
        self.redcorr = self.datacube['REDCORR'].data
        noise = np.ma.masked_values(self.datacube['IVAR'].data, 0).filled(1e-6)
        self.noise = np.sqrt(1/noise)
        # self.masked_flux = self.flux*np.logical_not(self.mask)
        if load_fit:
            self.model = np.ma.array(
                    self.datacube['MODEL'].data, 
                    mask=self.bm.flagged(self.datacube['MASK'].data, 
                                         ['FITIGNORED',]))
            self.stellarcontinuum = np.ma.array(    \
                    self.datacube['MODEL'].data     \
                    - self.datacube['EMLINE'].data  \
                    - self.datacube['EMLINE_BASE'].data, \
                    mask=self.bm.flagged(self.datacube['MASK'].data, \
                                         ['FITIGNORED',])) 
            self.emlines = np.ma.array( 
                    self.datacube['EMLINE'].data, 
                    mask=self.bm.flagged(self.datacube['EMLINE_MASK'].data, 
                                         ['ELIGNORED',]))
            self.emline_base = np.ma.array(
                    self.datacube['EMLINE_BASE'].data, 
                    mask=self.bm.flagged(self.datacube['EMLINE_MASK'].data, 
                                         ['ELIGNORED',]))
            self.residual = self.flux - self.model

    def __getitem__(self, xy):
        return self.getSpaxel(x=xy[0], y=xy[1])

    def getSpaxel(self, x=None, y=None, ra=None, dec=None):
        # x, y = y, x
        spaxel = Spaxel(self.plateifu)
        spaxel.wave = self.wave
        spaxel.wave_rest = self.wave_rest
        if ra is not None or dec is not None:
            sp_coor = coordinates.SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
            radec_xy = wcs.utils.skycoord_to_pixel(sp_coor, self.wcs_cude)
            x = np.round(radec_xy[1]).astype(int)
            y = np.round(radec_xy[0]).astype(int)
        spaxel.flux = self.flux[:, x, y]
        spaxel.redcorr = self.redcorr
        spaxel.noise = self.noise[:, x, y]
        spaxel.load_fit = self.load_fit
        if self.load_fit:
            spaxel.model = self.model[:, x, y]
            spaxel.stellarcontinuum = self.stellarcontinuum[:, x, y]
            spaxel.emlines = self.emlines[:, x, y]
            spaxel.emline_base = self.emline_base[:, x, y]
            spaxel.residual = self.residual[:, x, y]
        return spaxel


    # def __getitem__(self, item):
        # spaxel = Spaxel(self.plateifu)
        # spaxel.wave = self.wave
        # spaxel.wave_rest = self.wave_rest
        # spaxel.flux = self.flux[item]
        # spaxel.redcorr = self.redcorr
        # spaxel.noise = self.noise[item]
        # spaxel.load_fit = self.load_fit
        # if self.load_fit:
            # spaxel.model = self.model[item]
            # spaxel.stellarcontinuum = self.stellarcontinuum[item]
            # spaxel.emlines = self.emlines[item]
            # spaxel.emline_base = self.emline_base[item]
            # spaxel.residual = self.residual[item]
        # return spaxel
    
    def stackArea(self, areaSelect):
        spec_stacked = np.sum(np.moveaxis(self.flux, 0, -1)[areaSelect], 
                              axis=0)
        noise_stacked = np.sqrt(np.sum(np.moveaxis(
                                    self.noise**2, 0, -1)[areaSelect], axis=0))
        self.flux = spec_stacked
        self.noise = noise_stacked
        return spec_stacked

    def stack(self, binMatrix):
        if binMatrix is None:
            raise ValueError("bin number matrix must be given!")
        for binnum in np.unique(binMatrix):
            # flux = np.moveaxis(self.flux, 0, -1)
            # stellarcontinuum = np.moveaxis(self.stellarcontinuum, 0, -1)
            # noise = np.moveaxis(self.noise, 0, -1)
            binRegion = (binMatrix == binnum)
            self.flux[..., binRegion] = (np.sum(self.flux[..., binRegion], axis=1) 
                                            / np.sum(binRegion))[..., np.newaxis]
            self.stellarcontinuum[..., binRegion] = (np.sum(self.stellarcontinuum[..., binRegion], axis=1) / np.sum(binRegion))[..., np.newaxis]
            self.noise[..., binRegion] = (np.sqrt(np.sum(self.noise[..., binRegion]**2, axis=1)) / np.sum(binRegion))[..., np.newaxis]
            # self.flux = np.moveaxis(flux, -1, 0)
            # self.stellarcontinuum = np.moveaxis(stellarcontinuum, -1, 0)
            # self.noise = np.moveaxis(noise, -1, 0)

    def stack2(self, binmatrix, vmap=None):
        '''stack selected spaxels into one spectrum

        Params:
            binmatrix: the selecting matrix, with selected spaxels with True others with False

        Return:
            class Spaxel
        '''
        if binmatrix.shape != self.flux.shape[1:]:
            raise ValueError('binmatrix shape error!')
        spaxel = Spaxel(self.plateifu)
        spaxel.wave = self.wave
        spaxel.z = self.z
        spaxel.wave_rest = self.wave_rest
        if vmap is None:
            spaxel.flux = np.sum(self.flux[..., binmatrix], axis=1)
            spaxel.stellarcontinuum = np.sum(self.stellarcontinuum[..., binmatrix], axis=1)
            spaxel.noise = np.sqrt(np.sum(self.noise[..., binmatrix]**2, axis=1))
        else: #need to correct the rotation before stacking
            if binmatrix.shape != vmap.shape:
                raise ValueError('binmatrix amd velocity map unmatch!')
            a, b = binmatrix.shape
            flux = 0
            noise2 = 0 #squar of noise
            for i in range(a):
                for j in range(b):
                    if binmatrix[i, j]:
                        wave = self.wave / np.exp(vmap[i, j] / 299792.485)
                        flux += np.interp(spaxel.wave, wave, self.flux[:, i, j]) 
                        noise2 += np.interp(spaxel.wave, wave, self.noise[:, i, j])**2
            spaxel.flux = flux
            spaxel.noise = np.sqrt(noise2)
        return spaxel

    def genBinMatrix(self, emline, targetSN, quiet=True, mask=None, **kwargs):
        flux = self.flux_map2(emline, mask=mask, **kwargs)
        select = ~mask & (flux[0]/flux[1] > 0.5)
        binMatrix = np.full(flux[0].shape, -1)
        map_flux = flux[0][select]
        map_err = flux[1][select]
        # print("SNR:", np.sum(flux[0][select])/np.sum(flux[1][select]))
        if np.sum(map_flux)/np.sum(map_err) <= targetSN:
            if not quiet:
                print("cannot reach global {}".format(targetSN), 'uniform Bin-id was used!')
            return binMatrix
        if np.min(map_flux/map_err) > targetSN:
            if not quiet:
                print("Global SNR reached, no bin needed!")
            return None
        binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
            self.maps.x[select], self.maps.y[select], map_flux, \
            map_err, targetSN, plot=0, quiet=1)
        binMatrix[select] = binNum
        return binMatrix

    def flux_map(self, filename, mask=None, tie_balmer=False, limit_doublets=False, quiet=True,
                 directory='./', auto=False, broad_balmer=800, broad_O3=600, 
                 binmatrix=None, fewer_balmer=False, save_cube=False):
        '''output the map imformation into a fits file

        Params:
            filename: the output file name
        '''
        if mask is None:
            mask = self.maps.mask
        npix, naxis1, naxis2 = self.flux.shape
        lam_range_gal = np.array([np.min(self.wave), np.max(self.wave)])/(1+self.z)
        # get the template imformation
        velscale = 299792.458 * np.log(self.wave[1]/self.wave[0])
        FWHM_gal = 2.7
        miles = miles_util.miles('./ppxf_old/miles_models/Mun1.30*.fits', velscale, FWHM_gal)
        # get the fitted emission line
        gas_templates, gas_names, line_wave = ppxf_util.emission_lines(
                    miles.log_lam_temp, lam_range_gal, FWHM_gal, 
                    tie_balmer=tie_balmer, limit_doublets=limit_doublets, 
                    broad_balmer=True, broad_O3=True, quiet=True)
        multimap = dict(zip(gas_names, np.zeros((len(gas_names), 4, naxis1, naxis2))))
        datacube = np.zeros((*self.flux.shape, 3)) # 3 for balmer, forbidden, broad
        if binmatrix:
            self.stack(binmatrix)

        for x in range(naxis1):
            for y in range(naxis2):
                if mask[x, y]:
                    continue
                if quiet==2 or quiet==0:
                    print('coordinate:', [x, y])
                sp = self[x, y]
                if auto:
                    try:
                        try: 
                            # fitting with strong AGN
                            pp1 = sp.ppxf_fit(mode='emline', quiet=quiet, broad_balmer=800)
                            pp2 = sp.ppxf_fit(mode='emline', quiet=quiet, broad_balmer=800, 
                                              broad_O3=600)
                        except:
                            # fitting with weak AGN
                            if not quiet:
                                print('Change emline with fewer emission lines!')
                            pp1 = sp.ppxf_fit(mode='emline', quiet=quiet, fewer_lines=True)
                            pp2 = sp.ppxf_fit(mode='emline', quiet=quiet, fewer_lines=True, 
                                              broad_O3=600)
                    except KeyboardInterrupt:
                        sys.exit()
                    except: # for failed fitting
                        if not quiet:
                            print("Fitting failed!")
                        continue

                    F = (pp1.chi2_orig - pp2.chi2_orig)*(pp2.vars_num - pp2.params_num ) \
                            / (pp2.params_num - pp1.params_num) / pp2.chi2_orig
                    p_value = 1 - stats.f.cdf(F, pp2.params_num - pp1.params_num, 
                                                 pp2.vars_num - pp2.params_num)
                    if (p_value < 0.05) and ((pp2.chi2 - 1) < (pp1.chi2 - 1)) \
                            and np.any(pp2.gas_flux[-2:]/pp2.gas_flux_error[-2:] > 3):
                        pp = pp2
                        if not quiet:
                            print('Prefer broad [O III]')
                    else:
                        pp = pp1
                    if not quiet:
                        print('p_value:', p_value, 'fit1 chi2:', pp1.chi2, 
                              'fit2 chi2:', pp2.chi2)
                else:
                    pp = sp.ppxf_fit(mode='emline', quiet=quiet, broad_balmer=broad_balmer, 
                                     broad_O3=broad_O3)
                dwave = np.roll(pp.lam, -1) - pp.lam
                dwave[-1] = dwave[-2] # fix the bad point introduced by roll
                flux = dwave @ pp.matrix * pp.weights * pp.flux_scale
                flux_err = dwave @ pp.matrix \
                           * capfit.cov_err(pp.matrix / pp.noise[:, None])[1] \
                           * pp.flux_scale
                
                gas_flux = dict(zip(pp.gas_names, flux))
                gas_flux_err = dict(zip(pp.gas_names, flux_err))
                v, sigma = np.transpose(np.array(pp.sol)[pp.component.tolist()])
                rel_v = dict(zip(pp.gas_names, v - 299792.485 * np.log(1+pp.z)))
                sigma = dict(zip(pp.gas_names, sigma))
                for name in pp.gas_names:
                    multimap[name][:, x, y] = gas_flux[name], gas_flux_err[name],\
                                              rel_v[name], sigma[name]
                if save_cube:
                    for comp_n in [0, 1, 2]: # balmer, forbidden, and broad lines
                        wave_mask = (self.wave >= pp.lam[0]) & (self.wave <= pp.lam[-1])
                        comp_select = np.where(pp.component == comp_n)
                        cube_comp = pp.matrix[:, comp_select] @ (
                                    pp.weights[comp_select] * pp.flux_scale)
                        datacube[wave_mask, x, y, comp_n] = cube_comp[:, 0]
        hdr = fits.Header()
        hdr['AUTHER'] = 'cjhang'
        hdr['COMMENT'] = "Fitting emission lines with broad lines"
        primary_hdu = fits.PrimaryHDU(header=hdr)
        hdu_list = [primary_hdu]
        for name in gas_names:
            hdu_list.append(fits.ImageHDU(multimap[name], name=name))
        hdus = fits.HDUList(hdu_list)
        hdus.writeto('{}/{}.fits'.format(directory, filename), overwrite=True)

        ### saving the maps and cubes seperately
        # hdr = fits.Header()
        # hdr['AUTHER'] = 'cjhang'
        # hdr['COMMENT'] = "Fitting emission lines with broad lines"
        # primary_hdu = fits.PrimaryHDU(header=hdr)
        
        # # saveing the map data
        # hdu_maps = [primary_hdu]
        # for name in gas_names:
            # hdu_maps.append(fits.ImageHDU(multimap[name], name=name))

        # hdu_maps = fits.HDUList(hdu_maps)
        # hdu_maps.writeto('{}/{}-maps.fits'.format('./data', filename), overwrite=True)

        # # saving the fitted datacube data
        # # datacubes = np.zeros(np.unique(pp.component), *self.flux.shape)
        if save_cube:
            hdu_cubes_wave = fits.ImageHDU(self.wave, name='wave')
            hdu_cubes_flux = fits.ImageHDU(self.flux.data - self.stellarcontinuum.data, name='emline')
            hdu_cubes_fits = fits.ImageHDU(datacube, name='fits')
            hdu_cubes = [primary_hdu, hdu_cubes_wave, hdu_cubes_flux, hdu_cubes_fits]
            hdu_cubes = fits.HDUList(hdu_cubes)
            hdu_cubes.writeto('{}/{}-cubes.fits'.format(directory, filename), overwrite=True)

        # return multimap
        
    def flux_map2(self, names, mask=None, window=20, fix_invalid=True, **kwargs):
        m = Maps(self.plateifu)
        if isinstance(names, str):
            names = [names]
        #agn, *others = m.bptregion(snr=3)
        lines_num = len(names)
        a, b = m.naxis1, m.naxis2
        if mask is None:
            mask = m.mask
        flux_map = np.full((lines_num*8, a, b), 0) + 1e-8
        for i in range(a):
            for j in range(b):
                if not mask[i, j]:
                    sp = self[i, j]
                    #names = ['Hb-4862', 'OIII-5008', 'Ha-6564']
                    # print(i,j)
                    for k in range(lines_num):
                        emline = names[k]
                        fline1 = sp.fitline(emline, window=window, quiet=True, **kwargs) 
                        fline2 = sp.fitline(emline, window=window, broad_component=True, 
                                            quiet=True, **kwargs)
                        F = (fline1.chi2 - fline2.chi2)*(fline2.var_num - fline2.para_num ) \
                                / (fline2.para_num - fline1.para_num) / fline2.chi2
                        p_value = 1 - stats.f.cdf(F, fline2.para_num - fline1.para_num, 
                                                     fline2.var_num - fline2.para_num)
                        if (fline2.bn_ratio > 0.1) and (fline2.broad_snr > 3) \
                            and (fline2.broad_amp > 3*np.abs(fline2.noise).mean()) \
                            and (fline2.narrow_snr > 3) and (p_value < 0.05) \
                            and ((fline2.simple_chi2 - 1) < (fline1.simple_chi2 - 1)):
                            flux_map[8*k:8*k+8, i, j] = fline2.narrow_flux, \
                                                        fline2.narrow_flux_err, \
                                                        fline2.narrow_sigma, \
                                                        fline2.narrow_v, \
                                                        fline2.broad_flux, \
                                                        fline2.broad_flux_err, \
                                                        fline2.broad_sigma, \
                                                        fline2.broad_v
                        else:
                            # without fit the broad component
                            flux_map[8*k:8*k+8, i, j] = fline1.narrow_flux, \
                                                        fline1.narrow_flux_err, \
                                                        fline1.narrow_sigma, \
                                                        fline1.narrow_v, \
                                                        0, 0, 0, 0
        #flux_map[0][flux_map[0]/(flux_map[1]+1e-8) < 0.1] = -1 # remove the unreal pixels
        if fix_invalid:
            return utils.fixmap(flux_map)
        else:
            return flux_map
    
    def check_fit(self, emline, window=20, fit_mode='auto', filename=None):
        """
        check the spectrum fitting of a given line
        fit_mode: 
          - narrow, only narrow component
          - broad, add the broad component
          - auto, auto detect
        """
        m = Maps(self.plateifu)
        agn, *others = m.bptregion(snr=3)
        a, b = agn.shape
        if filename:
            pdf = mpdf.PdfPages(filename)
        else:
            today = date.today().isoformat()
            pdf = mpdf.PdfPages('results/checkfit-{}.pdf'.format(today))
        for i in range(a):
            for j in range(b):
                if agn[i, j]:
                    sp = self[:, i, j]
                    if fit_mode == 'narrow':
                        fig = sp.fitline(emline, window=window, quiet=True, 
                                            plot=True, return_fig=True) 
                    elif fit_mode == 'broad':
                        fig = sp.fitline(emline, window=window, broad_component=True, quiet=True,
                                            plot=True, return_fig=True)
                        
                    elif fit_mode == 'auto':
                        fline1 = sp.fitline(emline, window, quiet=True) 
                        fline2 = sp.fitline(emline, window, broad_component=True, quiet=True)
                        F = (fline1.chi2 - fline2.chi2)*(fline2.var_num - fline2.para_num ) \
                                / (fline2.para_num - fline1.para_num) / fline2.chi2
                        p_value = 1 - stats.f.cdf(F, fline2.para_num - fline1.para_num, 
                                                     fline2.var_num - fline2.para_num)
                        if (fline2.bn_ratio > 0.1) and (fline2.broad_snr > 3) \
                            and (fline2.narrow_snr > 3) and (p_value < 0.05) \
                            and ((fline2.simple_chi2 - 1) < (fline1.simple_chi2 - 1)):
                            fig = sp.fitline(emline, window, broad_component=True, quiet=True,
                                                plot=True, return_fig=True)
                        else:
                            fig = sp.fitline(emline, window, quiet=True, 
                                                plot=True, return_fig=True) 
                    plt.suptitle('[{}, {}]'.format(i,j))
                    pdf.savefig(fig)
                    plt.close()
    
        pdf.close()
        print("Done! Plot into ", 'results/checkfit-{}.pdf'.format(today))

