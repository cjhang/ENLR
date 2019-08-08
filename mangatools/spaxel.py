import re
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.table import Table
from astropy import wcs, coordinates

from .ppxf import ppxf, ppxf_util, miles_util, capfit

from .manga import MaNGA, package_path
from .maps import Maps

class Spaxel(MaNGA):
    # used for spaxels in datacube
    def __init__(self, plateifu):
        super().__init__(plateifu)
        self.fitted = False
        self.dx = None
        self.instru_sigma = None
        self.wave = []
        self.wave_rest = []
        self.flux = []
        self.noise = []
        # model data
        self.model = []
        self.stellarcontinuum = []
        self.emlines = []
        self.emline_base = []
        self.residual = []
        self.redcorr = None
        #self.sigma = self.psf/2.355

    def plot(self, waveRange=None, restWave=False, showFlux=True, showModel=False,
             ax=None, showEmlines=False, showContinuum=False, 
             showResidual=False, showImage=True):
        """show the spectrum, model and emission line
        """
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if restWave:
            wave = self.wave_rest
            ax.set_xlabel("wavelength (rest frame)")
        else:
            wave = self.wave
            ax.set_xlabel("wavelength")
        if showFlux:
            ax.step(wave, self.flux, label='flux', color='k', lw=0.5)
        if self.fitted or self.load_fit:
            if showModel:
                ax.plot(self.wave, self.model, label='model', color='r', lw=1)
            if showEmlines:
                ax.plot(self.wave, self.emlines, label='emlines', color='b', lw=1)
            if showContinuum:
                ax.plot(self.wave, self.stellarcontinuum, label='stellar continuum',
                        color='g', lw=1)
            if showResidual:
                ax.step(self.wave, self.residual-0.5, label='residual', 
                        where='mid', color='0.5', lw=0.5)
        if waveRange:
            wave_window = (wave > waveRange[0]) & (wave < waveRange[1]) 
            ax.set_xlim(waveRange)
            ax.set_ylim(self.flux[wave_window].min(), 
                        self.flux[wave_window].max())
        ax.set_ylabel('Flux/($10^{-17}ergs^{-1}cm^{-2}\AA^{-1}$)')
        ax.legend()

        if showImage:
            return fig
            #plt.show(fig)

    def ppxf_fit(self, tie_balmer=False, limit_doublets=False, mode='population', 
                 quiet=False, broad_balmer=None, broad_O3=None, fewer_lines=False):
        """fitting the spectrum using ppxf
        """
        # Only use the wavelength range in common between galaxy and stellar library.
        wave_mask = (self.wave > 3585*(1+self.z)) & (self.wave < 7340*(1+self.z))
        self.flux_scale = np.ma.median(self.flux[wave_mask])
        flux = self.flux[wave_mask]/self.flux_scale
        if self.redcorr is not None:
            flux = flux * self.redcorr[wave_mask]
        wave = self.wave[wave_mask]
        noise = self.noise[wave_mask]/self.flux_scale
        if not np.all((noise > 0) & np.isfinite(noise)) and not quiet:
            print('noise:', noise)
            print('flux_scale:', self.flux_scale)

        # noise = np.full_like(flux, 0.0166)
        c = 299792.485
        # define the dispersion of one pixel
        velscale = c*np.log(wave[1]/wave[0]) 
        FWHM_gal = self.psf
        # load the model
        miles = miles_util.miles(package_path + '/ppxf/miles_models/Mun1.30*.fits', velscale, FWHM_gal)
        reg_dim = miles.templates.shape[1:]
        stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)
        ## before fitting
        dv = c*(miles.log_lam_temp[0] - np.log(wave[0]))
        vel = c*np.log(1+self.z)
        start = [vel, 100.]
        if mode == 'population':
            regul_err = 0.013 # Desired regularization error
            # lam_range_gal = np.array([np.min(wave), np.max(wave)])
            lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1+self.z)
            gas_templates, gas_names, line_wave = ppxf_util.emission_lines(
                    miles.log_lam_temp, lam_range_gal, FWHM_gal,
                    tie_balmer=tie_balmer, limit_doublets=limit_doublets, 
                    broad_balmer=broad_balmer, broad_O3=broad_O3,
                    fewer_lines=fewer_lines)
            templates = np.column_stack([stars_templates, gas_templates])
            n_temps = stars_templates.shape[1]
            # Balmer lines start with contain letter
            n_balmer = np.sum([re.match('^[a-zA-Z]+', a) is not None for a in gas_names]) 
            # forbidden lines contain "["
            # n_forbidden = np.sum(["[" in a for a in gas_names])
            n_forbidden = np.sum([re.match('^\[[a-zA-Z]+', a) is not None 
                                    for a in gas_names])
            component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden 
            component_n = 3
            # for stars + balmer + forbidden lines
            moments = [4, 2, 2]
            start = [start, start, start]
            inf = 1e8
            bounds = [[[-inf, inf], [0, 500], [-inf, inf], [-inf, inf]], 
                      [[-inf, inf], [0, 500]], 
                      [[-inf, inf], [0, 500]]]
            if broad_balmer is not None:
                moments.append(2)
                start.append([vel, broad_balmer])
                # balmer up to 10000
                bounds.append([[-inf, inf], [broad_balmer, 10000]])
                # broad lines contain "{*}"
                n_balmer_broad = np.sum([re.match('^_[a-zA-Z]+', a) is not None 
                                            for a in gas_names]) 
                component = component + [component_n]*n_balmer_broad
                component_n = component_n + 1
            
            if broad_O3 is not None:
                moments.append(2)
                start.append([vel, broad_O3])
                bounds.append([[-inf, inf], [broad_O3, 2000]])
                n_forbidden_broad = np.sum([re.match('^_\[[a-zA-Z]+', a) is not None 
                                              for a in gas_names])
                component = component + [component_n]*n_forbidden_broad
            # print("moments:", moments)
            gas_component = np.array(component) > 0
            gas_reddening = 0 if tie_balmer else None
            # start fitting
            mask = ~flux.mask
            flux = flux.filled(1e-8)
            # noise = np.abs(noise.filled(np.ma.median(noise)))
            pp = ppxf.ppxf(templates, flux, noise, velscale, start,
                      plot=False, moments=moments, degree=12, mdegree=0, 
                      vsyst=dv, lam=wave, clean=False, regul=1/regul_err, 
                      reg_dim=reg_dim, mask=mask, component=component, 
                      quiet=quiet, gas_component=gas_component, 
                      gas_names=gas_names, gas_reddening=gas_reddening)
       
        elif mode == "kinematics":
            templates = stars_templates
            lamRange_temp = (np.exp(miles.log_lam_temp[0]), 
                            np.exp(miles.log_lam_temp[-1]))
            flux = np.ma.array(flux)
            bad_mask = np.where((flux.mask == True) 
                                & (flux > 10*np.ma.median(flux)))[0]
            goodpixels = ppxf_util.determine_goodpixels(
                            np.log(wave), lamRange_temp, self.z, 
                            broad_balmer=broad_balmer,
                            broad_O3=broad_O3)
            #goodpixels = np.flatnonzero(np.ones_like(flux))
            goodpixels = utils.deleteElement(goodpixels, bad_mask)
            #import pdb; pdb.set_trace();
            flux = flux.filled(0)
            pp = ppxf.ppxf(templates, flux, noise, velscale, start,
                    goodpixels=goodpixels, plot=False, moments=4, 
                    degree=12, vsyst=dv, clean=False, lam=wave, quiet=quiet)
            self.fitted = True
            self.stellarcontinuum = pp.bestfit * self.flux_scale
            self.wave_fit = pp.lam
            self.window_fit = ((self.wave <= self.wave_fit.max()) 
                                & (self.wave >= self.wave_fit.min()))
            self.emlines = self.flux[self.window_fit] - self.stellarcontinuum
        elif mode == 'emline':
            flux = (self.flux - self.stellarcontinuum) / self.flux_scale * self.redcorr
            flux = flux[wave_mask]
            # lam_range_gal = np.array([np.min(wave), np.max(wave)])
            lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1+self.z)
            gas_templates, gas_names, line_wave = ppxf_util.emission_lines(
                    miles.log_lam_temp, lam_range_gal, FWHM_gal, fewer_lines=fewer_lines, 
                    tie_balmer=tie_balmer, limit_doublets=limit_doublets, 
                    broad_balmer=broad_balmer, broad_O3=broad_O3, quiet=quiet)
            # Balmer lines start with contain letter
            n_balmer = np.sum([re.match('^[a-zA-Z]+', a) is not None for a in gas_names]) 
            # forbidden lines contain "["
            # n_forbidden = np.sum(["[" in a for a in gas_names])
            n_forbidden = np.sum([re.match('^\[[a-zA-Z]+', a) is not None 
                                    for a in gas_names])
            component = [0]*n_balmer + [1]*n_forbidden 
            # print('len component:', len(component), 'leb gas template:', len(gas_templates))
            component_n = 2
            # for stars + balmer + forbidden lines
            moments = [2, 2]
            start = [start, start]
            inf = 1e6
            bounds = [[[-inf, inf], [0, 600]], 
                      [[-inf, inf], [0, 600]]]
            if broad_balmer is not None:
                moments.append(2)
                start.append([vel, broad_balmer])
                # balmer up to 10000
                bounds.append([[-inf, inf], [broad_balmer, 10000]])
                # broad lines contain "{*}"
                n_balmer_broad = np.sum([re.match('^_[a-zA-Z]+', a) is not None 
                                            for a in gas_names]) 
                component = component + [2]*n_balmer_broad
                component_n = component_n + 1
            
            if broad_O3 is not None:
                moments.append(2)
                start.append([vel, broad_O3])
                bounds.append([[-inf, inf], [broad_O3, 2000]])
                n_forbidden_broad = np.sum([re.match('^_\[[a-zA-Z]+', a) is not None 
                                              for a in gas_names])
                component = component + [component_n]*n_forbidden_broad
            # print("moments:", moments)
            # print("component:", component)
            gas_component = np.array(component) >= 0
            gas_reddening = 0 if tie_balmer else None
            # start fitting
            mask = ~flux.mask
            # flux = flux.filled(1e-8)
            # print("flux shape:", flux.shape, 'noise shape:', noise.shape)
            # noise = np.ma.masked_invalid(noise) + 1e-8
            # flux = np.ma.masked_invalid(flux).filled(0)
            # import pdb; pdb.set_trace()
            flux = flux.filled(0)
            pp = ppxf.ppxf(gas_templates, flux, noise, velscale, start, bounds=bounds,
                      plot=False, moments=moments, degree=-1, mdegree=0, 
                      vsyst=dv, lam=wave, clean=False, mask=mask, component=component, 
                      quiet=quiet, gas_component=gas_component, 
                      gas_names=gas_names, gas_reddening=gas_reddening)
        # wrap relivant imformation into ppxf 
        pp.flux_scale = self.flux_scale
        pp.z = self.z
        pp.mode = mode
        # pp.var_num = len(pp.lam)
        # pp.para_num = len(np.concatenate(start))
        pp.chi2_orig = pp.chi2 * (pp.vars_num - pp.params_num)
        if False:#len(pp.gas_component) > 0:
            dwave = np.roll(pp.lam, -1) - pp.lam
            dwave[-1] = dwave[-2] # fix the bad point introduced by roll
            # dwave = np.ones_like(pp.lam)
            # TODO: combined flux_scale into weights
            pp.gas_flux = dwave @ pp.matrix * pp.weights * pp.flux_scale
            pp.gas_flux_err = (dwave @ pp.matrix 
                    * capfit.cov_err(pp.matrix / pp.noise[:, None])[1] * pp.flux_scale)
            pp.gas_lines = dict(zip(pp.gas_names, pp.gas_flux))
            pp.gas_lines_err = dict(zip(pp.gas_names, pp.gas_flux_error))

        return pp
        #return PPXFwrap(pp, self.flux_scale, self.z, mode=mode) 

    def autofit(self, mode='emline', quiet=1, broad_balmer=800, broad_O3=600, **kwargs):
        """auto choose whether adding the broad components
        """
        try:
            try: 
                # fitting with strong AGN
                pp1 = sp.ppxf_fit(mode='emline', quiet=quiet, broad_balmer=800)
                pp2 = sp.ppxf_fit(mode='emline', quiet=quiet, broad_balmer=800, 
                                  broad_O3=600)
            except: # if fitting failed
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
            raise ValueError("Fitting failed")

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
        return pp

    @staticmethod
    def gaussian(t, p):
        """
        Multiple Gaussian funtion
        # p[0]: amplitude
        # p[1]: mean
        # p[2]: sigma
        """
        y = np.zeros_like(t)
        for i in range(0, len(p), 3):
            amp = p[i]
            mean = p[i+1]
            sigma = p[i+2]
            y = y + amp * np.exp(-(t-mean)**2/(2*sigma**2))
        return y

    @staticmethod
    def cost_fun(p, t, y, y_err):
        return (Spaxel.gaussian(t, p) - y) / np.abs(y_err)

    @staticmethod
    def gaussian_list(t, p):
        y_list = []
        for i in range(0, len(p), 3):
            amp = p[i]
            mean = p[i+1]
            sigma = p[i+2]
            y_list.append(amp * np.exp(-(t-mean)**2/(2*sigma**2)))
        return np.array(y_list)

    @staticmethod
    def find_lines(emline, window):
        # return all the window contain all the lines
        # known_names = np.array(['Hb', '[OIII]', '[NII]_1', 'Ha', '[NII]_2'])
        known_lines = np.array([4862.69, 5008.24, 6549.84, 6564.61, 6585.23])
        near_lines = known_lines[(known_lines < emline + 2*window) \
                     & (known_lines > emline - 2*window)]
        wave_range = np.array([np.min(near_lines) - window, np.max(near_lines) + window])
        return near_lines, wave_range

    @staticmethod
    def gen_err(t, dt, p, p_err, time=1000):
        # generate gaussian profile error
        new_p = np.zeros((p.size, time))
        for i in range(len(p)):
            new_p[i,:] = p[i] + np.random.randn(time) * p_err[0]
        # expand to have same dimension with t, to do broadcast calculation
        new_p2 = np.repeat(new_p[:, :, np.newaxis], t.size, axis=2)
        ret_err = np.std(np.sum(Spaxel.gaussian(t, new_p2) * dt, axis=1))
        return ret_err

