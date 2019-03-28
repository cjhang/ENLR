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
        if self.fitted:
            if showModel:
                ax.plot(self.wave_fit, self.model, label='model', color='r', lw=1)
            if showEmlines:
                ax.plot(self.wave_fit, self.emlines, label='emlines', color='b', lw=1)
            if showContinuum:
                ax.plot(self.wave_fit, self.stellarcontinuum, label='stellar continuum',
                        color='g', lw=1)
            if showResidual:
                ax.step(self.wavefit, self.residual-0.5, label='residual', 
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
        wave_mask = (self.wave > 3540*(1+self.z)) & (self.wave < 7400*(1+self.z))
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
        FWHM_gal = 2.7
        # load the model
        miles = miles_util.miles(package_path + '/ppxf/miles_models/Mun1.30*.fits', velscale, FWHM_gal)
        reg_dim = miles.templates.shape[1:]
        stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)
        ## before fitting
        dv = c*(miles.log_lam_temp[0] - np.log(wave[0]))
        vel = c*np.log(1+self.z)
        start = [vel, 180.]
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

    def fitline(self, name, window=20, velocity=150, broad_component=False, 
                broad_velocity=500, plot=False, quiet=False, return_fig=False, 
                debug=False):
        """
        Fiting the emission line using gaussian profile
        This procedure try to self-determine how many components in this lines
        """
        # fitting preparation
        names = ['Hb-4862', 'OIII-5008', 'Ha-6564']
        wavelengths = [4862.69, 5008.24, 6564.61]
        emline = dict(zip(names, wavelengths))[name]

        # check double lines 
        if name == 'Ha-6564':
            return self.fitHa(broad_component=broad_component, 
                              broad_velocity=broad_velocity, plot=plot, 
                              quiet=quiet, return_fig=return_fig, debug=debug)
        lines, wave_range = Spaxel.find_lines(emline, window)
        wave_window = (self.wave_rest > wave_range[0]) \
                      & (self.wave_rest < wave_range[1])
        wave = self.wave_rest[wave_window] - emline
        dwave = (np.roll(self.wave_rest, -1) - self.wave_rest)[wave_window]
        flux = ((self.flux - self.stellarcontinuum) * self.redcorr)[wave_window]
        # fixed those region with all masked
        flux = np.ma.masked_invalid(flux).filled(0) 
        flux_err = self.noise[wave_window]
        scale = 1 #np.median(self.flux)
        # spectrum parameters
        dx = np.mean(wave[1:] - wave[:-1])
        logdx = (np.log(self.wave[-1] + emline) \
                - np.log(wave[0] + emline))/(wave.size - 1)
        #instrument sigma in Angstrom
        instru_sigma = self.psf / 2.355 / logdx / emline 
        sigma = velocity * emline / (self.c * dx)
        broad_sigma = broad_velocity * emline / (self.c * dx)
        line_window = (wave > -window/3) & (wave < window/3)
        amp = np.abs(np.max(flux[line_window])) / scale
        if broad_component:
            #assert broad_velocity
            # two components of the line
            x0 = np.array([amp, 0, sigma, amp/10, 0, broad_sigma]) 
            bound_low = np.array([0, -window/3, instru_sigma, 0, -window, 
                                 broad_sigma])
            bound_up = np.array([np.inf, window/3, broad_sigma, np.inf, 
                                 window, 3000])
        else:
            x0 = np.array([amp, 0, sigma]) # two components of the line
            bound_low = np.array([0, -window/3, instru_sigma])
            bound_up = np.array([np.inf, window/3, broad_sigma])
        relative_lines = lines - emline
        near_lines = relative_lines[np.abs(relative_lines) > np.min(np.abs(
                relative_lines))]
        for line in near_lines:
            line_window = (wave > line - window/3) & (wave < line + window/3)
            amp = np.max(flux[line_window]) / scale
            x0 = np.append(x0, [amp, line, sigma])
            bound_low = np.append(bound_low, [0, line - window/3, instru_sigma])
            bound_up = np.append(bound_up, [np.inf, line + window/3, broad_sigma])
        def cost_fun(p, t, y, y_err):
            return (Spaxel.gaussian(t, p) - y) / np.abs(y_err)
        if debug:
            bounds = list(zip(bound_low, bound_up))
            for i in range(len(x0)):
                if (x0[i] < bounds[i][0]) or (x0[i]>bounds[i][1]):
                    print("----------Initial value issue!------------")
                    print('x0:', x0)
                    print('bound_low:', bound_low)
                    print('bound_up:', bound_up)

        res_lsq = optimize.least_squares(Spaxel.cost_fun, x0, 
                                         args=(wave, flux/scale, flux_err/scale), \
                                         bounds=[bound_low, bound_up], jac='3-point')
        # derived parameters
        perror = capfit.cov_err(res_lsq.jac)[1]
        chi2 = np.sum((flux - scale * Spaxel.gaussian(wave, res_lsq.x))**2 / flux_err**2)
        flux_fit = scale * Spaxel.gaussian_list(wave, res_lsq.x)
        # return values
        class Line(): pass
        ret_line = Line()
        # res_lsq.x [amp, mean, sigma]
        ret_line.noise = flux_err
        ret_line.narrow_amp = res_lsq.x[0]
        ret_line.narrow_v = res_lsq.x[1] * (self.c * dx) / emline
        ret_line.narrow_sigma = res_lsq.x[2] * (self.c * dx) / emline
        # theorectical integrate S = amp * \sqrt{2*pi*sigma^2}
        # print('flux broad sum up:', np.sum(flux_fit[1] * dwave))
        ret_line.narrow_flux = res_lsq.x[0]*np.sqrt(2*np.pi*res_lsq.x[2]**2)
        ret_line.narrow_flux_err = scale * Spaxel.gen_err(wave, dwave, res_lsq.x[:3], perror[:3])
        ret_line.narrow_snr = ret_line.narrow_flux / (ret_line.narrow_flux_err + ESP)
        if broad_component:
            ret_line.broad_amp = res_lsq.x[3]
            ret_line.broad_v = res_lsq.x[4] * (self.c * dx) / emline
            ret_line.broad_sigma = res_lsq.x[5] * (self.c * dx) / emline
            ret_line.broad_flux = res_lsq.x[3]*np.sqrt(2*np.pi*res_lsq.x[5]**2)
            ret_line.broad_flux_err = scale * Spaxel.gen_err(wave, dwave, res_lsq.x[3:6], perror[3:6])
            ret_line.broad_snr = ret_line.broad_flux / (ret_line.broad_flux_err + ESP)
            ret_line.bn_ratio = ret_line.broad_flux / (ret_line.narrow_flux + ESP)
            # maximum flux ratio
            ret_line.bn_amp_ratio = np.max(flux_fit[1]) / (np.max(flux_fit[0])+ESP) 
        ret_line.chi2 = chi2
        ret_line.var_num = wave.size
        ret_line.para_num = res_lsq.x.size
        ret_line.simple_chi2 = chi2 / (wave.size - res_lsq.x.size)
        if not quiet:
            print("==========================================")
            for i in range(0, len(x0), 3):
                print("Component:", i//3)
                print("Initial value: {:.2f} {:.2f} {:.2f}".format(*x0[i:i+3]))
                print("Fitting bounds:")
                print("  lower: {low[0]:.2f} {low[1]:.2f} {low[2]:.2f}".format(low=bound_low[i:i+3]))
                print("  upper: {up[0]:.2f} {up[1]:.2f} {up[2]:.2f}".format(up=bound_up[i:i+3]))
                print("Parameters: {:.2f} {:.2f} {:.2f}".format(*res_lsq.x[i:i+3]))
                print("Errors: {:.2f} {:.2f} {:.2f}".format(*perror[i:i+3]))
                print("-----------------------------------------")
            # calculate the Chi^2
            print("Chi^2:", ret_line.simple_chi2)
            print("narrow component SNR:", ret_line.narrow_snr)
            if broad_component:
                print("broad-narrow amplitude ratio:", ret_line.bn_amp_ratio)
                print("broad-narrow flux ratio:", ret_line.bn_ratio)
                print("broad component SNR:", ret_line.broad_snr)
            print("==========================================")
        if plot:
            fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
            plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2, wspace=0.1, hspace=0.0)
            ax1.step(wave + emline, flux, label='data', where='mid')
            if len(flux_fit) > 1:
                for i in range(len(flux_fit)):
                    ax1.plot(wave + emline, flux_fit[i], label='component {}'.format(i))
            ax1.plot(wave+emline, scale * Spaxel.gaussian(wave, res_lsq.x), label='fit')
            ax1.set_xticklabels([])
            ax1.legend()
            ax1.set_ylabel('flux')
            ax2.plot(wave+emline, flux - scale * Spaxel.gaussian(wave, res_lsq.x), label='residual')
            ax2.set_xlabel('wavelength')
        if return_fig:
            return fig
        else:
            return ret_line
    
    def fitHa(self, window=50, velocity=150, broad_component=False, broad_velocity=500, 
                plot=False, quiet=False, return_fig=False, debug=False):
        """
        Fiting the emission line using gaussian profile
        This procedure try to self-determine how many components in this lines
        """
        # fitting preparation
        emline = 6564.61
        wave_range = emline + np.array([-window, window])
        wave_window = (self.wave_rest > wave_range[0]) \
                      & (self.wave_rest < wave_range[1])
        wave = self.wave_rest[wave_window] - emline
        dwave = (np.roll(self.wave_rest, -1) - self.wave_rest)[wave_window]
        flux = ((self.flux - self.stellarcontinuum) * self.redcorr)[wave_window]
        flux = np.ma.masked_invalid(flux).filled(0)
        flux_err = self.noise[wave_window]
        scale = 1 #np.median(self.flux)
        # spectrum parameters
        dx = np.mean(wave[1:] - wave[:-1])
        logdx = (np.log(self.wave[-1] + emline) - np.log(wave[0] + emline))/(wave.size - 1)
        instru_sigma = self.psf / 2.355 / logdx / emline #instrument sigma in Angstrom
        sigma = velocity * emline / (self.c * dx)
        rot_vel = 400 * emline / (self.c * dx)
        broad_sigma = broad_velocity * emline / (self.c * dx)
        Ha_line_window = (wave > -window/3) & (wave < window/3)
        Ha_amp = np.abs(np.max(flux[Ha_line_window])) / scale
        N2_line_window = (wave > 20.62-window/3) & (wave < 20.62+window/3)
        N2_amp = np.abs(np.max(flux[N2_line_window])) / scale
        lines = np.array([6549.84, 6564.61, 6585.23]) - emline
        x0 = [Ha_amp, 0, sigma, N2_amp, sigma]
        bound_low = np.array([0, -rot_vel, instru_sigma, 0, instru_sigma])
        bound_up = np.array([np.inf, rot_vel, broad_sigma, np.inf, broad_sigma])
        # model the lines
        def fit_narrow(t, p):
            y_fit = p[0] * np.exp(-(t-lines[1]-p[1])**2/(2*p[2]**2)) + \
                    p[3]/3 * np.exp(-(t-lines[0]-p[1])**2/(2*p[4]**2)) + \
                    p[3] * np.exp(-(t-lines[2]-p[1])**2/(2*p[4]**2))
            return y_fit
        def fit_broad(t, p):
            y_fit = p[0] * np.exp(-(t-lines[1]-p[1])**2/(2*p[2]**2)) + \
                    p[3]/3 * np.exp(-(t-lines[0]-p[1])**2/(2*p[4]**2)) + \
                    p[3] * np.exp(-(t-lines[2]-p[1])**2/(2*p[4]**2)) + \
                    p[-3] * np.exp(-(t-p[-2])**2/(2*p[-1]**2))
            return y_fit
        if not broad_component:
            fit_fun = fit_narrow
        else:
            x0 = np.append(x0, [Ha_amp/10, 0, 2*broad_sigma]) # two components of the line
            bound_low = np.append(bound_low, [0, -window, broad_sigma])
            bound_up = np.append(bound_up, [np.inf, window, 3000])
            fit_fun = fit_broad
        def cost_fun(p, t, y, y_err):
            y_fit = fit_fun(t, p)
            return (y_fit - y) / np.abs(y_err)
        if debug:
            bounds = list(zip(bound_low, bound_up))
            for i in range(len(x0)):
                if (x0[i] < bounds[i][0]) or (x0[i]>bounds[i][1]):
                    print("----------Initial value issue!------------")
                    print('x0:', x0)
                    print('bound_low:', bound_low)
                    print('bound_up:', bound_up)

        res_lsq = optimize.least_squares(cost_fun, x0, args=(wave, flux/scale, 
                                         flux_err/scale), 
                                         bounds=[bound_low, bound_up], jac='3-point')
        # derived parameters
        perror = capfit.cov_err(res_lsq.jac)[1]
        chi2 = np.sum((flux - scale * fit_fun(wave, res_lsq.x))**2 / flux_err**2)
        # return values
        class Line(): pass
        ret_line = Line()
        # res_lsq.x [amp, mean, sigma]
        ret_line.noise = flux_err
        ret_line.narrow_amp = res_lsq.x[0]
        ret_line.narrow_v = res_lsq.x[1] * (self.c * dx) / emline
        ret_line.narrow_sigma = res_lsq.x[2] * (self.c * dx) / emline
        ret_line.narrow_fit = Spaxel.gaussian(wave, res_lsq.x[:3])
        # ret_line.narrow_flux = np.sum(ret_line.narrow_fit * dwave)
        ret_line.narrow_flux = res_lsq.x[0] * np.sqrt(2*np.pi*res_lsq.x[2]**2)
        ret_line.narrow_flux_err = scale * Spaxel.gen_err(wave, dwave, res_lsq.x[:3], perror[:3])
        ret_line.narrow_snr = ret_line.narrow_flux / (ret_line.narrow_flux_err + ESP)
        if broad_component:
            ret_line.broad_amp = res_lsq.x[-3]
            ret_line.broad_v = res_lsq.x[-2] * (self.c * dx) / emline
            ret_line.broad_sigma = res_lsq.x[-1] * (self.c * dx) / emline
            ret_line.broad_fit = Spaxel.gaussian(wave, res_lsq.x[-3:])
            # ret_line.broad_flux = np.sum(ret_line.broad_fit * dwave)
            ret_line.broad_flux = res_lsq.x[-3] * np.sqrt(2*np.pi*res_lsq.x[-1]**2)
            ret_line.broad_flux_err = scale * Spaxel.gen_err(wave, dwave, res_lsq.x[-3:], perror[-3:])
            ret_line.broad_snr = ret_line.broad_flux / (ret_line.broad_flux_err + ESP)
            ret_line.bn_ratio = ret_line.broad_flux / (ret_line.narrow_flux + ESP)
            ret_line.bn_amp_ratio = np.max(ret_line.broad_fit) / (np.max(ret_line.narrow_fit) + ESP)
        ret_line.chi2 = chi2
        ret_line.var_num = wave.size
        ret_line.para_num = res_lsq.x.size
        ret_line.simple_chi2 = chi2 / (wave.size - res_lsq.x.size)
        if not quiet:
            print("==========================================")
            print("Narrow component:")
            print("Initial value:", x0)
            print("Fitting bounds:")
            print("  lower: {}".format(bound_low))
            print("  upper: {}".format(bound_up))
            print("Parameters: {}".format(*res_lsq.x))
            print("Errors: {}".format(*perror))
            print("-----------------------------------------")
            # calculate the Chi^2
            print("Chi^2:", ret_line.simple_chi2)
            print("narrow component SNR:", ret_line.narrow_snr)
            if broad_component:
                print("broad-narrow amplitude ratio:", ret_line.bn_amp_ratio)
                print("broad-narrow flux ratio:", ret_line.bn_ratio)
                print("broad component SNR:", ret_line.broad_snr)
            print("==========================================")
        if plot:
            fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
            plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2, wspace=0.1, hspace=0.0)
            ax1.step(wave + emline, flux, label='data', where='mid')
            if broad_component:
                ax1.plot(wave + emline, fit_narrow(wave, res_lsq.x[:-3]), label='narrow fit')
                ax1.plot(wave+emline, Spaxel.gaussian(wave, res_lsq.x[-3:]), label='broad fit')
                ax1.plot(wave+emline, fit_fun(wave, res_lsq.x), label='fit')
            else:
                ax1.plot(wave + emline, fit_narrow(wave, res_lsq.x), label='narrow fit')
            ax1.set_xticklabels([])
            ax1.legend()
            ax1.set_ylabel('flux')
            ax2.plot(wave+emline, flux - scale * fit_fun(wave, res_lsq.x), label='residual')
            ax2.set_xlabel('wavelength')
        if return_fig:
            return fig
        else:
            return ret_line

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

