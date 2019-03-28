import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.cosmology import LambdaCDM

from .config import SAS

class SDSS():
    """direct get the info and emline imformation of SDSS single fiber data
    """

    COSMO = LambdaCDM(H0=70, Om0=0.3, Ode0=0.6)

    ## for the SDSS dr7 data matched with manga
    SDSS_MaNGA = SAS + '/other/dr7-mpl6_v4.fits'
    try:
        with  fits.open(SDSS_MaNGA) as SDSS_DR7:
            DR7_INFO = SDSS_DR7[1].data
            DR7_LINE = SDSS_DR7[2].data
    except:
        raise ValueError("No manga data found!")

    def __init__(self, plateifu):
        self.plateifu = plateifu
        self.sdss = SDSS.DR7_INFO[SDSS.DR7_INFO['plateifu'] == self.plateifu]
        self.mangaid = self.sdss['mangaid'][0]
        self.z = self.sdss['z'][0]
        self.z_err = self.sdss['z_err'][0]
        self.d = SDSS.COSMO.comoving_transverse_distance(self.z)
        self.d_err = self.z_err/self.z*self.d
        if (self.sdss['angle_dist'][0] < 2/3600) and (
            self.sdss['kpc_dist'][0]*u.kpc < self.d_err):
            self.match = True
        else:
            self.match = False
    
    @property
    def info(self):
        return SDSS.DR7_INFO.columns.names

    def detail(self):
        # check the matching quality
        if (self.sdss['angle_dist'][0] > 1e-4) or (
            self.sdss['kpc_dist'][0]*u.kpc > self.d_err):
            print("plateifu:{} have bad matching quality".format(
                      self.plateifu))
        else:
            print('good matching')
    
    @property
    def line_info(self):
        em_lines = ['H_alpha', 'H_delta', 'OII_3726', 'OII_3729', 'NEIII_3869',
                    'H_gamma', 'OIII_4363', 'H_beta', 'OIII_4959', 'OIII_5007',
                    'HEI_5876', 'OI_6300', 'NII_6548', 'NII_6584', 'SII_6717',
                    'SII_6731', 'ARIII7135', ]
        return em_lines

    def line(self, name, redcorr=False):
        emlines = SDSS.DR7_LINE[SDSS.DR7_INFO['plateifu']==self.plateifu]
        line_obs = emlines[name+'_FLUX'][0] 
        line_obs_err = emlines[name+'_FLUX_ERR'][0]
        if redcorr:
            Ha = emlines['H_alpha_FLUX'][0]
            Hb = emlines['H_beta_FLUX'][0]
            line_wave = re.match(r"\w+_([\d, \w]+)", name)[1]
            if line_wave in ['alpha', 'beta', 'gamma', 'delta']:
                wave_dict = {'alpha':6562., 'beta':4861., 'gamma':4340, 
                             'delta':4101}
                line_wave = wave_dict[line_wave]
            if Ha < 1e-8 or Hb < 1e-8: # 1e-8 the error
                return np.nan, np.nan
            line_obs = utils.dust_corr(Ha, Hb, float(line_wave)/10000, 
                                          line_obs).data
        return np.array([line_obs, line_obs_err])
    
    def lineEW(self, name):
        emlines = SDSS.DR7_LINE[SDSS.DR7_INFO['plateifu']==self.plateifu]
        return -emlines[name+'_EQW'][0], emlines[name+'_EQW_ERR'][0]

    def bpt(self, strickMode=False, plot=False, fig=None, index=None, show_kewly06=True):
        #print(self.line('H_alpha')[0])
        #if self.line('H_alpha')[0] < 1e-8 or H_beta < 1e-8: # 1e-8 the error
            # print('Invalid line!')
        x0 = self.line('NII_6584')[0]/self.line('H_alpha')[0]
        y0 = self.line('OIII_5007')[0]/self.line('H_beta')[0]
        #z, z_err = self.lineEW('H_alpha') # it's color, EW(Ha)
        region_name = ['AGN', 'CP', 'SF', 'seyfert', 'liner']
        x = np.ma.array(np.log10(x0))
        y = np.ma.array(np.log10(y0))
        types = utils.bptregion(np.ma.array(x), np.ma.array(y))
        # print('{}: {}'.format(self.plateifu, types))
        if strickMode:
            for t in range(1,len(types)):
                if types[t]:
                    gtype = region_name[t]
                    break
            else:
                gtype = None
        else:
            for t in range(len(types)):
                if types[t]:
                    gtype = region_name[t]
                    break
            else:
                gtype = None
        if fig == None:
            fig = plt.figure()
        if index == None:
            index = 111
        ax = fig.add_subplot(index)
        if show_kewly06:
            # Ke01 line [Kewly et al (2006)]
            x1 = np.linspace(-2.0, 0.46, 50)
            y_ke01 = 0.61/(x1-0.47)+1.19
            # Ka03 line [Kauffmann et al.(2003)]
            x2 = np.linspace(-2.0, 0, 50)
            y_ka03 = 0.61/(x2-0.05)+1.3
            # Schawinski2007 imperical line
            x3 = np.linspace(-0.2, 0.46, 20)
            y_schawinski = 1.05*x3+0.45
            ax.plot(x1, y_ke01, 'r', zorder=1)
            ax.plot(x2, y_ka03, 'k--', zorder=2)
            ax.plot(x3, y_schawinski, 'k--', zorder=3)
            ax.set_xlim(-2.0, 1.5)
            ax.set_ylim(-1.5, 1.5)
            ax.set_xlabel(r'LOG([NII]/H$\alpha)$')
            ax.set_ylabel(r'LOG([OIII]/H$\beta)$')
        ax.plot(x, y, 'b*', markersize=10.0, label=gtype)
        if plot:
            plt.show(fig)
        else:
            return gtype
