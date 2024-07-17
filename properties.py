from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.profile import HaloDensityProfile
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties import PropertyCalculation, LivePropertyCalculation
import pynbody
import numpy as np


#Radial Momentum profile property for calculating in/outflow rates
def p_r_weighted(self, weight=None):
    p_r = np.zeros(self.nbins)
    cumind = []
    for i in range(self.nbins):
        subs = self.sim[self.binidi]
        if weight != None:
            weight = subs[weight]
        else:
            weight = 1
        p_r[i] = ((weight*subs['vr'].in_units('kpc yr**-1'))*subs['mass'].in_units('Msol')).sum()
    return p_r

# Mass Flux
@pynbody.analysis.profile.Profile.profile_property
def p_r(self):
    return p_r_weighted(self)

# Metal Mass Flux
@pynbody.analysis.profile.Profile.profile_property
def p_r_metals(self):
    return p_r_weighted(self, weight='metals')

# Oxygen Mass Flux
@pynbody.analysis.profile.Profile.profile_property
def p_r_Ox(self):
    return p_r_weighted(self, weight='OxMassFrac')

# Iron Mass Flux
@pynbody.analysis.profile.Profile.profile_property
def p_r_Fe(self):
    return p_r_weighted(self, weight='FeMassFrac')

class InflowOutflow(HaloDensityProfile):
    '''
    Inflow and Outflow rate calculations.
    '''
    names = 'mass_outflow_profile', 'mass_inflow_profile', 'metal_outflow_profile', 'metal_inflow_profile', 'Ox_outflow_profile', 'Ox_inflow_profile', 'Fe_outflow_profile', 'Fe_inflow_profile'

    def __init__(self, simulation):
        super().__init__(simulation)
        # Filters to select inflowing and outflowing particles.
        self.ifilt = pynbody.filt.LowPass('vr', 0)
        self.ofilt = pynbody.filt.HighPass('vr', 0)

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def plot_ylabel(self):
        return r"Mass Flux $(M_\odot yr^{-1})$"

    def _get_profile(self, halo, maxrad):
        delta = self.plot_xdelta()
        nbins = int(maxrad / delta)
        maxrad = delta * nbins
        halo['ones'] = np.ones(len(halo))

        oprof = pynbody.analysis.profile.Profile(halo.g[self.ofilt], type='lin', ndim=3,
                                           min=0, max=maxrad, nbins=nbins, weight='mass')
        iprof = pynbody.analysis.profile.Profile(halo.g[self.ifilt], type='lin', ndim=3,
                                           min=0, max=maxrad, nbins=nbins, weight='mass')
        mass_out = oprof['p_r']/delta
        mass_in = iprof['p_r']/delta
        metal_out = oprof['p_r_metals']/delta
        metal_in = iprof['p_r_metals']/delta
        Ox_out = oprof['p_r_Ox']/delta
        Ox_in = iprof['p_r_Ox']/delta
        Fe_out = oprof['p_r_Fe']/delta
        Fe_in = iprof['p_r_Fe']/delta
        return mass_out, mass_in, metal_out, metal_in, Ox_out, Ox_in, Fe_out, Fe_in
        
    @centred_calculation
    def calculate(self, halo, existing_properties):
        try:
            vcen = pynbody.analysis.halo.vel_center(halo,cen_size=maxrad, retcen=True)
        except:
            return None, None, None

        halo['vel'] -= vcen

        mass_out, mass_in, Z_out, Z_in, Ox_out, Ox_in, Fe_out, Fe_in = self._get_profile(halo, existing_properties["max_radius"])

        halo['vel'] += vcen
        return mass_out, mass_in, Z_out, Z_in, Ox_out, Ox_in, Fe_out, Fe_in

class MetalProfile(PynbodyPropertyCalculation):
    '''
    Mass weighted metalicity profiles
    '''
    names = 'star_metal_profile', 'gas_metal_profile', 'star_Fe_profile', 'gas_Fe_profile', 'star_Ox_profile', 'gas_Ox_profile'

    def __init__(self, simulation):
        super().__init__(simulation)
        self._xdelta = self.get_simulation_property("approx_resolution_kpc", 0.1)

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def plot_x0(self):
        return self.plot_xdelta()/2

    def plot_xdelta(self):
        return self._xdelta

    def plot_xlabel(self):
        return "r/kpc"

    def plot_ylabel(self):
        return r"Metal Mass Fraction"

    def _get_profile(self, halo, maxrad):
        delta = self.plot_xdelta()
        nbins = int(maxrad / delta)
        maxrad = delta * nbins

        metal_pro_s = None
        metal_pro_g = None
        fe_pro_g = None
        fe_pro_s = None
        ox_pro_g = None
        ox_pro_s = None

        if len(halo.g)>10:
            pro_g = pynbody.analysis.profile.Profile(halo.g, type='lin', ndim=3,
                                               min=0, max=maxrad, nbins=nbins, weight='mass')
            metal_pro_g = pro_g['metals']
            fe_pro_g = pro_g['FeMassFrac']
            ox_pro_g = pro_g['OxMassFrac']
        if len(halo.s)>10:
            pro_s = pynbody.analysis.profile.Profile(halo.s[pynbody.filt.HighPass('tform',0)], type='lin', ndim=3,
                                               min=0, max=maxrad, nbins=nbins, weight='mass')
            metal_pro_s = pro_s['metals']
            fe_pro_s = pro_s['FeMassFrac']
            ox_pro_s = pro_s['OxMassFrac']
        return metal_pro_g, metal_pro_s, fe_pro_g, fe_pro_s, ox_pro_g, ox_pro_s

    @centred_calculation
    def calculate(self, halo, existing_properties):
        metals_gas, metals_star, fe_gas, fe_star, ox_gas, ox_star = self._get_profile(halo, existing_properties["max_radius"])
        return metals_star, metals_gas, fe_star, fe_gas, ox_star, ox_gas


class ColdDenGasMetalProfile(PynbodyPropertyCalculation):
    '''
    Mass weighted metalicity profiles for cold gas only (T < 2e4 K)
    '''
    names = 'cold_gas_metal_profile', 'cold_gas_Fe_profile', 'cold_gas_Ox_profile'

    def __init__(self, simulation):
        super().__init__(simulation)
        self._xdelta = self.get_simulation_property("approx_resolution_kpc", 0.1)

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def plot_x0(self):
        return self.plot_xdelta()/2

    def plot_xdelta(self):
        return self._xdelta

    def plot_xlabel(self):
        return "r/kpc"

    def plot_ylabel(self):
        return r"Metal Mass Fraction"

    def _get_profile(self, halo, maxrad):
        delta = self.plot_xdelta()
        nbins = int(maxrad / delta)
        maxrad = delta * nbins

        metal_pro = None
        fe_pro = None
        ox_pro = None

        if len(halo)>10:
            pro = pynbody.analysis.profile.Profile(halo, type='lin', ndim=3,
                                               min=0, max=maxrad, nbins=nbins, weight='mass')
            metal_pro = pro['metals']
            fe_pro = pro['FeMassFrac']
            ox_pro = pro['OxMassFrac']

        return metal_pro, fe_pro, ox_pro

    @centred_calculation
    def calculate(self, halo, existing_properties):
        cold_metal_pro, cold_fe_pro, cold_ox_pro = self._get_profile(halo.g[pynbody.filt.LowPass('temp',2e4)],
                                           existing_properties['max_radius'])
        #dense_metal_pro, dense_fe_pro, dense_ox_pro = self._get_profile(halo.g[pynbody.filt.HighPass('rho', '0.1 m_p cm**-3')],
         #                                                            existing_properties['max_radius'])
        return cold_metal_pro, cold_fe_pro, cold_ox_pro#, dense_metal_pro, dense_fe_pro, dense_ox_pro


class MassEnclosedTemp(PynbodyPropertyCalculation):
    '''
    Enclosed mass for different temperature ranges (MIGHT NEED TO CHANGE THESE)
    '''
    @classmethod
    def name(self):
        return "cold_gas_mass_profile", "warm_gas_mass_profile", "hot_gas_mass_profile"
     
    def requires_property(self):
        return ["shrink_center", "max_radius"]
     
    def rstat(self, halo, maxrad,delta=0.1):
        nbins = int(maxrad / delta)
        maxrad = delta * (nbins + 1)
        if len(halo.g)==0:
            return None, None, None
        proCG = pynbody.analysis.profile.Profile(halo.g[pynbody.filt.LowPass("temp", 1.e5)],
                                                 type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
        proWG = pynbody.analysis.profile.Profile(halo.g[pynbody.filt.BandPass("temp", 1.e5, 1e6)],
                                                 type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
        proHG = pynbody.analysis.profile.Profile(halo.g[pynbody.filt.HighPass("temp", 1.e6)],
                                                 type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)

        return proCG['mass_enc'], proWG['mass_enc'], proHG['mass_enc']
     
    @centred_calculation
    def calculate(self,halo,properties):
        maxrad = properties['max_radius']
        delta = properties.get('delta',0.1)

        CG, WG, HG = self.rstat(halo,maxrad,delta)\

        return CG, WG, HG

    def plot_x0(cls):
        return 0.0

    def plot_xdelta(cls):
        return 0.1

class Radius(LivePropertyCalculation):
    '''
    calculates a radius relative to the cricial density (i.e. Radius(200) = R200 based on total mass profile)
    '''
    def __init__(self, simulation, n_crit=200):
        super(Radius, self).__init__(simulation)
        self._omegaM0 = 0.3086
        self._omegaL0 = 0.6914
        self._h0 = 0.67
        self._ncrit = n_crit
        
    @classmethod
    def name(cls):
        return 'radius'
    
    def requires_property(self):
        return ['tot_mass_profile']
    
    def live_calculate(self, halo,*args):
        ts = halo.timestep
        z = ts.redshift
        a = 1.0 / (1.0 + z)
        H_z = pynbody.analysis.cosmology._a_dot(a, self._h0, self._omegaM0, self._omegaL0) / a
        H_z = pynbody.units.Unit("100 km s^-1 Mpc^-1") * H_z
        rho_crit = (3 * H_z ** 2) / (8 * np.pi * pynbody.units.G)
        if halo['tot_mass_profile'].max() == 0:
            return None
        else:
            rho_mean = halo['tot_mass_profile']/(4./3. * np.pi * ((np.arange(len(halo['tot_mass_profile']))+1)*0.1)**3)
            return np.where(rho_mean>rho_crit.in_units('Msol kpc**-3')*self._ncrit)[0][-1]*0.1

class StellarProfileFaceOn(PynbodyPropertyCalculation):
    '''
    calculate surface brightness
    '''
    @classmethod
    def name(self):
        return "v_surface_brightness", "b_surface_brightness", "i_surface_brightness"

    def plot_x0(cls):
        return 0.05

    @classmethod
    def plot_xdelta(cls):
        return 0.1

    @classmethod
    def plot_xlabel(cls):
        return "R/kpc"

    @classmethod
    def plot_ylog(cls):
        return False

    @classmethod
    def plot_xlog(cls):
        return False

    @staticmethod
    def plot_ylabel():
        return "v mags/arcsec$^2$", "b mags/arcsec$^2$", "i mags/arcsec$^2$"

    def requires_property(self):
        return ['shrink_center', 'max_radius']

    def calculate(self, halo, existing_properties):
        with pynbody.analysis.angmom.faceon(halo, cen = existing_properties['shrink_center']):
            nbins = int(existing_properties['max_radius']/0.1)
            ps = pynbody.analysis.profile.Profile(halo.s, type='lin', ndim=2, min=0, max=existing_properties['max_radius'], nbins=nbins)
            vals = [ps['sb,'+x] for x in ('v','b','i')]
        return vals
    
    
class StellarProfileDiagnosis(LivePropertyCalculation):
    '''
    calculate sersic profile fit properties based on existing surface brightness profiles
    '''
    def __init__(self, simulation, band, sats=0, sblimit=28, type='hlr', smooth=0):
        '''
        :param band: which band of calculated sb profile to use (band_surface_brightness)
        :param sats: if 1, search for nearby halos and limit the maximum radius to fit over accordingly
        :param type: hlr = limit max radius to 5 times half light radius. sb = limit based on sblimit
        :param sblimit: lowest surface brightness to consider (will cut off once sb reaches below this limit, only if type='sb')
        :param smooth: smooth over this number of bins when fitting. Useful when resolution limit is courser than spatial bin size
        '''
        super(StellarProfileDiagnosis, self).__init__(simulation)
        self.band = band
        self.type=type
        self.sblimit=sblimit
        self.sats=sats
        self.smooth=smooth
        if smooth==0:
            self.smooth=1

    names="half_light","sersic_m0", "sersic_n", "sersic_r0"

    @classmethod
    def requires_particle_data(self):
        return False

    def requires_property(self):
        return ["v_surface_brightness", "b_surface_brightness", "i_surface_brightness", 'max_radius', 'shrink_center']
    
    def sersic_surface_brightness(r, mueff, reff, n):
        # I(R) = m0 + 2.5*b_n/ln(10) ( (r/reff)^(1/n)-1 )
        #b_n taken based on solution to gamma function (Capaccioli 1989) for n < 10.
        return mueff + 2.5*(0.868*n-0.142)*((r/reff)**(1./n) - 1)

    def fit_sersic(r, surface_brightness, return_cov=False):
        s0_guess = np.mean(surface_brightness[:3])
        #s0_range = [s0_guess-2, s0_guess+2]
        #r0_range = [0.1,r[-1]]
        #n_range = [0.1, 6.0]
        s0_range = [10,40]
        n_range = [0.5,16.5]
        r0_range=[0,100]

        r0_guess = min(r[int(len(r)/2)], 50)
        n_guess = 1.0

        sigma = 10**(0.6*(surface_brightness-20))/r
        sigma = None

        popt, pcov = scipy.optimize.curve_fit(sersic_surface_brightness,r,surface_brightness,
                                          bounds=np.array((s0_range, r0_range, n_range)).T,
                                          sigma=sigma,
                                          p0=(s0_guess, r0_guess, n_guess))

        if return_cov:
            return popt, pcov
        else:
            return popt

    def live_calculate(self, halo, *args):
        r0 = 0.05
        delta_r = self.get_simulation_property("aprox_resolution_kpc", 0.1)
        if self.band+"_surface_brightness" in halo.keys():
            surface_brightness = halo[self.band+"_surface_brightness"]
        else:
            try:
                surface_brightness = halo.calculate(self.band+"_surface_brightness")
            except NoResultsError:
                return None, None, None, None
        if len(surface_brightness)<self.smooth*2:
            return None, None, None, None
        if surface_brightness.max()==0:
            return None, None, None, None
        r = np.arange(len(surface_brightness))*delta_r + 0.05
        nbins = len(surface_brightness)
        if self.smooth > 1:
            surface_brightness_new = np.nanmean(np.pad(surface_brightness[self.smooth:].astype(float),
                                                (0,3-surface_brightness[self.smooth:].size%self.smooth),mode='constant',
                                                constant_values=np.NaN).reshape(-1,self.smooth),axis=1)
            surface_brightness = np.insert(surface_brightness_new,0,np.mean(surface_brightness[:self.smooth]))
            r = np.arange(len(surface_brightness))*(0.1*self.smooth) + 0.1*self.smooth/2.
            r[0] = r0


        flux_density = 10**(surface_brightness/-2.5)
        flux_density[flux_density!=flux_density]=0
        #nbins = len(surface_brightness)
        #r = np.arange(0,halo['max_radius']+delta_r,delta_r)
        cumu_flux_density = (r * flux_density).cumsum()
        cumu_flux_density/=cumu_flux_density[-1]

        try:
            half_light_i = np.where(cumu_flux_density>0.5)[0][0]
            half_light = r0+delta_r * half_light_i * self.smooth
        except:
            half_light = None
        if self.type=='hlr':
            if half_light is None:
                return None, None, None, None
            maxrad = half_light*5

        if self.type=='sb':
            maxradi = np.where(surface_brightness>self.sblimit)[0]
            if len(maxradi) == 0:
                maxrad = halo['max_radius']
            else:
                maxrad = r[maxradi[0]]

        if self.sats:
            cen, mvir = halo.timestep.gather_property('shrink_center', 'finder_mass')
            darr = halo['shrink_center'] - cen
            D = np.sqrt(np.sum(darr**2,axis=1))
            dmin = D[(D>0)].min()
            if dmin < maxrad:
                maxrad = dmin

        usefit = np.where(r<=maxrad)[0]
        r_fit = r[usefit]
        sb_fit = surface_brightness[usefit]

        mask_not_nan = sb_fit==sb_fit
        r_fit = r_fit[mask_not_nan]
        sb_fit = sb_fit[mask_not_nan]
        if len(r_fit)<2:
            return None, None, None, None
        else:
            try:
                m0, r0, n = fit_sersic(r_fit, sb_fit)
            except:
                m0 = None
                n = None
                r0 = None
            return half_light, m0, n, r0

class VelDispersionProfile(PynbodyPropertyCalculation):
    '''
    dispersion of stars, gas and dm profiles
    '''
    @classmethod
    def name(self):
        return "vrdisp_stars", "vrdisp_gas", "vrdisp_dm", "vdisp_stars_3d", "vdisp_gas_3d", "vdisp_dm_3d"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def rstat(self, halo, maxrad,delta=0.1):
        nbins = int(maxrad / delta)
        maxrad = delta * (nbins + 1)
        sigg = None
        sigs = None
        sigdm = None
        sigs3D = None
        sigg3D = None
        sigdm3D = None
        if len(halo.g)>5:
            proG = pynbody.analysis.profile.Profile(halo.g, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
            sigg =  proG['vr_disp']
            sigg3D = np.sqrt(proG['vx_disp']**2 + proG['vy_disp']**2 + proG['vz_disp']**2)
        if len(halo.s)>5:
            proS = pynbody.analysis.profile.Profile(halo.s, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
            sigs = proS['vr_disp']
            sigs3D = np.sqrt(proS['vx_disp'] ** 2 + proS['vy_disp'] ** 2 + proS['vz_disp'] ** 2)
        if len(halo.dm)>5:
            proDM = pynbody.analysis.profile.Profile(halo.dm, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
            sigdm = proDM['vr_disp']
            sigdm3D = np.sqrt(proDM['vx_disp'] ** 2 + proDM['vy_disp'] ** 2 + proDM['vz_disp'] ** 2)
        return sigs, sigg, sigdm, sigs3D, sigg3D, sigdm3D

    @centred_calculation
    def calculate(self,halo,properties):
        maxrad = properties['max_radius']
        delta = properties.get('delta',0.1)
        try:
            vcen = pynbody.analysis.halo.vel_center(halo,cen_size='5 kpc', retcen=True)
        except:
            return None, None, None

        halo['vel'] -= vcen

        sigS, sigG, sigDM, sigS3D, sigG3D, sigDM3D = self.rstat(halo,maxrad,delta)

        halo['vel'] += vcen

        return sigS, sigG, sigDM, sigS3D, sigG3D, sigDM3D

    def plot_x0(cls):
        return 0.0

    def plot_xdelta(cls):
        return 0.1

#Define two profile properties used below
@pynbody.analysis.profile.Profile.profile_property
def vr_disp_encl(self):
    vrdisp = np.zeros(self.nbins)
    cumind = []
    for i in range(self.nbins):
        cumind = np.append(cumind,self.binind[i])
        cumind = np.sort(cumind).astype(np.int)
        subs = self.sim[cumind]
        # use = np.where(subs.g['temp'] > temp_cut)[0]
        vr = subs['vr'].view(np.ndarray)
        weight = subs[self._weight_by].view(np.ndarray)
        mean_sq = ((vr*weight).sum()/weight.sum())**2
        sq_mean = (vr**2*weight).sum()/weight.sum()
        vrdisp[i] = np.sqrt(sq_mean - mean_sq)
    return vrdisp

@pynbody.analysis.profile.Profile.profile_property
def v_disp_tot_encl(self):
    v3d_disp = np.zeros(self.nbins)
    cumind = []
    for i in range(self.nbins):
        cumind = np.append(cumind,self.binind[i])
        cumind = np.sort(cumind).astype(np.int)
        subs = self.sim[cumind]
        # use = np.where(subs.g['temp'] > temp_cut)[0]
        vx = subs['vx'].view(np.ndarray)
        vy = subs['vy'].view(np.ndarray)
        vz = subs['vz'].view(np.ndarray)
        weight = subs[self._weight_by].view(np.ndarray)
        mean_sq_x = ((vx*weight).sum()/weight.sum())**2
        sq_mean_x = (vx**2*weight).sum()/weight.sum()
        sigx2 = sq_mean_x - mean_sq_x
        mean_sq_y = ((vy * weight).sum() / weight.sum()) ** 2
        sq_mean_y = (vy** 2 * weight).sum() / weight.sum()
        sigy2 = sq_mean_y - mean_sq_y
        mean_sq_z = ((vz * weight).sum() / weight.sum()) ** 2
        sq_mean_z = (vz ** 2 * weight).sum() / weight.sum()
        sigz2 = sq_mean_z - mean_sq_z
        v3d_disp[i] = np.sqrt(sigx2 + sigy2 + sigz2)
    return v3d_disp

class VelDispersionProfileEncl(PynbodyPropertyCalculation):
    '''
    enclosed velosity dispersion profiles
    '''
    @classmethod
    def name(self):
        return "vrdisp_encl_stars", "vrdisp_encl_gas", "vrdisp_encl_dm", "vrdisp_encl_stars_3d", "vdisp_encl_gas_3d", "vrdisp_encl_dm_3d"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def rstat(self, halo, maxrad,delta=0.1):
        nbins = int(maxrad / delta)
        maxrad = delta * (nbins + 1)
        sigG = None
        sigS = None
        sigDM = None
        sigG3d = None
        sigS3d = None
        sigDM3d = None
        proG = pynbody.analysis.profile.Profile(halo.g, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
        proS = pynbody.analysis.profile.Profile(halo.s, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)
        proDM = pynbody.analysis.profile.Profile(halo.dm, type='lin', ndim=3, min=0, max=maxrad, nbins=nbins)

        if len(halo.g)>5:
            sigG = proG['vr_disp_encl']
            sigG3d = proG['v_disp_tot_encl']
        if len(halo.s)>5:
            sigS = proS['vr_disp_encl']
            sigS3d = proG['v_disp_tot_encl']
        if len(halo.dm)>5:
            sigDM = proDM['vr_disp_encl']
            sigDM3d = proG['v_disp_tot_encl']

        return sigS, sigG, sigDM, sigS3d, sigG3d, sigDM3d

    @centred_calculation
    def calculate(self,halo,properties):
        maxrad = properties['max_radius']
        delta = properties.get('delta',0.1)
        try:
            vcen = pynbody.analysis.halo.vel_center(halo,cen_size='5 kpc', retcen=True)
        except:
            return None, None, None

        halo['vel'] -= vcen

        sigS, sigG, sigDM = self.rstat(halo,maxrad,delta)

        halo['vel'] += vcen

        return sigS, sigG, sigDM
