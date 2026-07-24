import functools
from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.profile import HaloDensityProfile
from tangos.properties.pynbody.BH import BHAccHistogram
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties import LivePropertyCalculation
from tangos.live_calculation import NoResultsError
import pynbody
import numpy as np
import scipy.optimize

# Metal species mass fraction arrays, in the order their profiles are stored.
METAL_ARRAYS = ('metals', 'FeMassFrac', 'OxMassFrac')

# Photometric bands, in the order their surface brightness profiles are stored.
SB_BANDS = ('u', 'g', 'r', 'i', 'z', 'U', 'V', 'J')

# Specific angular momentum profile keys, in the order their profiles are stored.
ANGMOM_KEYS = ('jtot', 'j_phi', 'j_theta')

# Radial flux profile properties and the per particle array each weights the
# radial momentum by, in the order their in/outflow properties are stored.
FLUX_PROFILES = {'mass': ('p_r', None),
                 'metal': ('p_r_metals', 'metals'),
                 'Ox': ('p_r_Ox', 'OxMassFrac'),
                 'Fe': ('p_r_Fe', 'FeMassFrac')}

# Temperature bands for the enclosed gas mass profiles, in stored order.
TEMPERATURE_BANDS = {'cold': pynbody.filt.LowPass("temp", 1.e5),
                     'warm': pynbody.filt.BandPass("temp", 1.e5, 1.e6),
                     'hot': pynbody.filt.HighPass("temp", 1.e6)}

# A species is only profiled if it has more than this many particles.
MIN_PARTICLES = 10

# The radius vel_center averages the bulk halo velocity over. Every velocity
# centred calculation below uses this same value.
VEL_CENTER_SIZE = '5 kpc'


# Replace NaN values with a substitute value in some list of arrays. Species
# that were not profiled at all are passed through as None.
def nan_remover(arrs, sub=0):
    return [arr if arr is None else np.nan_to_num(arr, nan=sub) for arr in arrs]


def null_result(prop):
    '''A result of the right shape for a calculation that cannot be made'''
    if isinstance(prop.names, str):
        return None
    return (None,) * len(prop.names)


def lin_profile(region, maxrad, nbins, ndim=3, **kwargs):
    '''A linearly binned pynbody profile of region, running from 0 out to maxrad'''
    return pynbody.analysis.profile.Profile(region, type='lin', ndim=ndim,
                                            min=0, max=maxrad, nbins=nbins, **kwargs)


def family_profile(family, maxrad, nbins, select=None, min_particles=MIN_PARTICLES, **kwargs):
    '''
    As lin_profile, but None for a species with too few particles to say anything
    useful about. The particle count is checked explicitly because catching
    pynbody's exceptions for empty families is expensive. Pass select to profile
    only a subset of the family while still counting the family as a whole.
    '''
    if len(family) <= min_particles:
        return None
    return lin_profile(family if select is None else select, maxrad, nbins, **kwargs)


def profile_values(profiles, keys):
    '''
    Pull each key from each profile, ordered key by key (the order in which the
    properties below name their arrays), substituting None for species that were
    not profiled.
    '''
    return [None if pro is None else pro[key] for key in keys for pro in profiles]


def load_metal_arrays(halo):
    '''Force the metal arrays of both gas and stars to load before profiling'''
    for family in (halo.g, halo.s):
        for array in METAL_ARRAYS:
            family[array]


def hot_phase_key(halo):
    '''
    The name of the hot phase gas mass array, or None if this snapshot has no
    two phase gas. Which of the two spellings is used depends on the loader.
    '''
    for key in ('massHot', 'MassHot'):
        try:
            halo.g[key]
        except Exception:
            continue
        return key
    return None


def two_phase_keys(halo):
    '''
    The names of the (hot, cold) phase gas mass arrays, or (None, None) if this
    snapshot has no two phase gas. Only the hot phase mass is stored, so the
    cold phase array is derived here as mass - massHot.
    '''
    hot = hot_phase_key(halo)
    if hot is None:
        return None, None
    cold = hot[:-3] + 'Cold'
    halo.g[cold] = halo.g['mass'] - halo.g[hot]
    return hot, cold


def velocity_centred(family=None):
    '''
    Decorator for calculate() methods that must work in the halo's rest frame.

    The bulk velocity found by pynbody's vel_center over the central
    VEL_CENTER_SIZE is subtracted before the calculation and added back
    afterwards, since the snapshot is shared between property calculations.
    family restricts the shift to a single particle family, so that velocities
    the calculation does not use are not needlessly rounded. If the bulk
    velocity cannot be found, the whole property group is skipped.
    '''
    def decorator(calculate):
        @functools.wraps(calculate)
        def calculate_in_rest_frame(self, halo, existing_properties):
            try:
                vcen = pynbody.analysis.halo.vel_center(halo, return_cen=True,
                                                        cen_size=VEL_CENTER_SIZE)
            except Exception:
                return null_result(self)
            region = halo if family is None else getattr(halo, family)
            region['vel'] -= vcen
            try:
                return calculate(self, halo, existing_properties)
            finally:
                region['vel'] += vcen
        return calculate_in_rest_frame
    return decorator


class BinnedProfileProperty(HaloDensityProfile):
    '''
    Shared machinery for the profile properties below. names is None so that
    tangos does not register this helper as a calculation of the density
    profiles it would otherwise inherit.
    '''
    names = None

    def binning(self, maxrad):
        '''
        Bin width, number of bins, and maxrad rounded down to a whole number of
        bins. Every profile here shares this binning so that the stored arrays
        line up bin for bin.
        '''
        delta = self.plot_xdelta()
        nbins = int(maxrad / delta)
        return delta, nbins, delta * nbins


#Radial Momentum profile property for calculating in/outflow rates
@pynbody.analysis.profile.Profile.profile_property
def p_r_weighted(self, weight=None):
    vr_array = self.sim['vr'].view(np.ndarray)
    mass_array = self.sim['mass'].view(np.ndarray)
    if weight is not None:
        weight_array = self.sim[weight].view(np.ndarray)
    else:
        weight_array = np.ones(vr_array.size)
    product = (weight_array*vr_array*mass_array)
    result = np.array([product[self.binind[i]].sum() for i in range(self.nbins)])
    result = result.view(pynbody.array.SimArray)
    result.units = self.sim['vr'].units*self.sim['mass'].units
    result.sim = self.sim
    return result.in_units('Msol kpc yr**-1')


def _register_flux_profile(name, weight):
    '''Register a profile property giving the radial momentum weighted by weight'''
    def flux(self):
        return p_r_weighted(self, weight=weight)
    flux.__name__ = name
    return pynbody.analysis.profile.Profile.profile_property(flux)


# Mass flux (p_r), and the metal, oxygen and iron mass fluxes derived from it.
for _name, _weight in FLUX_PROFILES.values():
    _register_flux_profile(_name, _weight)


class InflowOutflow(BinnedProfileProperty):
    '''
    Inflow and Outflow rate calculations.
    '''
    names = 'mass_outflow_profile', 'mass_inflow_profile', 'metal_outflow_profile', 'metal_inflow_profile', 'Ox_outflow_profile', 'Ox_inflow_profile', 'Fe_outflow_profile', 'Fe_inflow_profile'

    # Filters to select inflowing and outflowing particles.
    ifilt = pynbody.filt.LowPass('vr', 0)
    ofilt = pynbody.filt.HighPass('vr', 0)

    def plot_ylabel(self):
        return r"Mass Flux $(M_\odot yr^{-1})$"

    def _get_profile(self, halo, maxrad):
        load_metal_arrays(halo)
        delta, nbins, maxrad = self.binning(maxrad)
        # Outflowing particles carry positive and inflowing particles negative
        # radial momentum; both are stored as positive rates.
        directions = ((1, lin_profile(halo.g[self.ofilt], maxrad, nbins)),
                      (-1, lin_profile(halo.g[self.ifilt], maxrad, nbins)))
        return tuple((sign * pro[key] / delta).view(np.ndarray)
                     for key, _ in FLUX_PROFILES.values()
                     for sign, pro in directions)

    @centred_calculation
    @velocity_centred(family='g')
    def calculate(self, halo, existing_properties):
        return self._get_profile(halo, existing_properties["max_radius"])


class AngMomProfile(BinnedProfileProperty):
    '''
    Mass weighted specific angular momentum magnitude (j) and direction (jphi/jtheta)
    '''
    names = 'j_gas_profile', 'j_star_profile', 'j_dm_profile',\
    'j_phi_gas_profile', 'j_phi_star_profile', 'j_phi_dm_profile',\
    'j_theta_gas_profile', 'j_theta_star_profile', 'j_theta_dm_profile'

    def plot_xlabel(self):
        return "r/kpc"

    def plot_ylabel(self):
        return (r"j_{tot} kpc km s$^{-1}$",)*3 + (r"$\phi_j$ [radians]",)*3 \
            + (r"$\theta_j$ [radians]",)*3

    def _get_profile(self, halo, maxrad):
        delta, nbins, maxrad = self.binning(maxrad)
        profiles = [family_profile(family, maxrad, nbins)
                    for family in (halo.g, halo.s, halo.dm)]
        return profile_values(profiles, ANGMOM_KEYS)

    @centred_calculation
    @velocity_centred()
    def calculate(self, halo, existing_properties):
        return self._get_profile(halo, existing_properties["max_radius"])


class MetalProfile(BinnedProfileProperty):
    '''
    Mass weighted metalicity profiles
    '''
    names = 'star_metal_profile', 'gas_metal_profile', 'star_Fe_profile', 'gas_Fe_profile', 'star_Ox_profile', 'gas_Ox_profile'

    def plot_xlabel(self):
        return "r/kpc"

    def plot_ylabel(self):
        return r"Metal Mass Fraction"

    def _get_profile(self, halo, maxrad):
        delta, nbins, maxrad = self.binning(maxrad)
        load_metal_arrays(halo)
        # Black holes are stored as star particles with a negative formation time
        formed = halo.s[pynbody.filt.HighPass('tform', 0)]
        profiles = (family_profile(halo.s, maxrad, nbins, select=formed, weight='mass'),
                    family_profile(halo.g, maxrad, nbins, weight='mass'))
        return profile_values(profiles, METAL_ARRAYS)

    @centred_calculation
    def calculate(self, halo, existing_properties):
        return tuple(nan_remover(self._get_profile(halo, existing_properties["max_radius"])))


class ColdDenGasMetalProfile(BinnedProfileProperty):
    '''
    Mass weighted metalicity profiles for cold gas only (T < 2e4 K)
    '''
    names = 'cold_gas_metal_profile', 'cold_gas_Fe_profile', 'cold_gas_Ox_profile'

    def plot_ylabel(self):
        return r"Metal Mass Fraction"

    def _get_profile(self, halo, maxrad):
        delta, nbins, maxrad = self.binning(maxrad)
        load_metal_arrays(halo)
        if len(halo) <= MIN_PARTICLES:
            return null_result(self)

        cold = pynbody.filt.LowPass('temp', 2e4)
        hot_key, cold_key = two_phase_keys(halo)
        if hot_key is None:
            pro = lin_profile(halo[cold], maxrad, nbins, weight='mass')
            return tuple(pro[array] for array in METAL_ARRAYS)

        # Two phase particles hold cold gas whatever their nominal temperature,
        # so they are profiled separately, weighted by their cold phase mass.
        twophase = pynbody.filt.HighPass(hot_key, 0)
        one_pro = lin_profile(halo[~twophase & cold], maxrad, nbins, weight='mass')
        cold_pro = lin_profile(halo[twophase], maxrad, nbins, weight=cold_key)
        return tuple(one_pro[array] + cold_pro[array] for array in METAL_ARRAYS)

    @centred_calculation
    def calculate(self, halo, existing_properties):
        return tuple(nan_remover(self._get_profile(halo.g,
                                                   existing_properties['max_radius'])))


class InstantaneousSFR(PynbodyPropertyCalculation):
    '''
    Calculate the instantaneous star formation rate based on the SFR of all eligible
    gas particles in the simulation.
    '''
    names = "instantaneous_SFR"
    def calculate(self, halo, existing_properties):
        # Load SF Parameters
        dCStar = self.get_simulation_property("dCStar", 0.05)
        dPhysDenMin = self.get_simulation_property("dPhysDenMin",
                                                   0.1)*pynbody.units.m_p/(pynbody.units.cm**3)
        dTempMax = self.get_simulation_property("dTempMax", 1.5e4)
        H2 = "COOLING_MOLECULARH" in self.get_simulation_property("macros", "")

        # Only select particles that are eligible to form stars
        # (rho > dPhysDenMin, temp < dTempMax, not two phase)
        dense_filter = pynbody.filt.HighPass('rho', dPhysDenMin)
        cold_filter = pynbody.filt.LowPass('temp', dTempMax)
        eligible = dense_filter & cold_filter
        hot_key = hot_phase_key(halo)
        if hot_key is not None:
            eligible &= ~pynbody.filt.HighPass(hot_key, 0)
        sfg = halo.g[eligible]

        # Calculate the instantaneous SFR for all gas particles.
        tdyn = 1.0/np.sqrt(4*np.pi*pynbody.units.G*sfg['rho'])
        sfr = dCStar*(sfg['mass']/tdyn).in_units('Msol yr**-1')
        if H2:
            return (sfr*sfg['H2']).sum()
        else:
            return sfr.sum()


class MassEnclosedTemp(BinnedProfileProperty):
    '''
    Enclosed mass for different temperature ranges (MIGHT NEED TO CHANGE THESE)
    '''
    names = "cold_gas_mass_profile", "warm_gas_mass_profile", "hot_gas_mass_profile"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def plot_ylabel(self):
        return "Gas Mass"

    def rstat(self, halo, maxrad):
        delta, nbins, maxrad = self.binning(maxrad)
        if len(halo.g)==0:
            return null_result(self)

        hot_key, cold_key = two_phase_keys(halo)
        if hot_key is None:
            return tuple(lin_profile(halo.g[band], maxrad, nbins)['mass_enc']
                         for band in TEMPERATURE_BANDS.values())

        # Two phase particles hold hot and cold gas simultaneously whatever their
        # nominal temperature, so they are profiled separately and each phase
        # mass added to the band it belongs in.
        twophase = pynbody.filt.HighPass(hot_key, 0)
        two_pro = lin_profile(halo.g[twophase], maxrad, nbins)
        phase_mass = {'cold': cold_key, 'hot': hot_key}
        enclosed = []
        for name, band in TEMPERATURE_BANDS.items():
            mass_enc = lin_profile(halo.g[~twophase & band], maxrad, nbins)['mass_enc']
            if name in phase_mass:
                mass_enc = mass_enc + pynbody.analysis.profile.weight_fn(two_pro,
                                                                         phase_mass[name]).cumsum()
            enclosed.append(mass_enc)
        return tuple(enclosed)

    @centred_calculation
    def calculate(self,halo,properties):
        return self.rstat(halo, properties['max_radius'])


class LiveRadius(LivePropertyCalculation):
    '''
    calculates a radius relative to the cricial density (i.e. Radius(200) = R200 based on total mass profile)
    '''
    names = 'radius'
    def __init__(self, simulation, n_crit=200):
        super(LiveRadius, self).__init__(simulation)
        self._ncrit = n_crit


    def requires_property(self):
        return ['star_mass_profile', 'gas_mass_profile','dm_mass_profile']

    def live_calculate(self, halo,*args):
        self._deltar = self.get_simulation_property("approx_resolution_kpc", 0.1)
        self._omegaM0 = self.get_simulation_property("dOmega0",0.3086)
        self._omegaL0 = self.get_simulation_property("dLambda", 0.6914)
        self._h0 = self.get_simulation_property("dHubble0", 2.89443)
        self._dunit = self.get_simulation_property("dKpcUnit", 25000)
        self._munit = self.get_simulation_property("dMsolUnit", 1.99101e+15)
        ts = halo.timestep
        z = ts.redshift
        a = 1.0 / (1.0 + z)
        denunit = self._munit/self._dunit**3
        velunit = 8.0285 * np.sqrt(6.6743e-8 * denunit) * self._dunit
        hubunit = 10. * velunit / self._dunit
        self._h0 *= hubunit

        H_z = pynbody.analysis.cosmology._a_dot(a, self._h0, self._omegaM0, self._omegaL0) / a
        H_z = pynbody.units.Unit("100 km s^-1 Mpc^-1") * H_z
        rho_crit = (3 * H_z ** 2) / (8 * np.pi * pynbody.units.G)

        tot_mass_profile = halo['dm_mass_profile']+halo['gas_mass_profile']+halo['star_mass_profile']
        rho_mean = tot_mass_profile/(4./3. * np.pi * ((np.arange(len(tot_mass_profile))+1)*self._deltar)**3)
        return (np.where(rho_mean>rho_crit.in_units('Msol kpc**-3')*self._ncrit)[0][-1]+1)*self._deltar


class StellarProfileFaceOn(BinnedProfileProperty):
    '''
    calculate surface brightness
    '''
    names = "u_surface_brightness", "g_surface_brightness", "r_surface_brightness", "i_surface_brightness", "z_surface_brightness", "U_surface_brightness", "V_surface_brightness", "J_surface_brightness"

    @classmethod
    def plot_xlog(cls):
        return False

    @classmethod
    def plot_ylog(cls):
        return False

    @staticmethod
    def plot_ylabel():
        return tuple(band+" mags/arcsec$^2$" for band in SB_BANDS)

    def requires_property(self):
        return ['shrink_center', 'max_radius']

    @centred_calculation
    @velocity_centred()
    def calculate(self, halo, existing_properties):
        delta, nbins, maxrad = self.binning(existing_properties['max_radius'])
        with pynbody.analysis.angmom.faceon(halo, already_centered=True):
            ps = lin_profile(halo.s, maxrad, nbins, ndim=2)
            return [ps['sb,'+band] for band in SB_BANDS]


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
        self.requires_particle_data = False
        if smooth==0:
            self.smooth=1

    names="half_light","sersic_m0", "sersic_n", "sersic_r0"

    def requires_property(self):
        return list(StellarProfileFaceOn.names) + ['max_radius', 'shrink_center']

    def sersic_surface_brightness(self,r, mueff, reff, n):
        # I(R) = m0 + 2.5*b_n/ln(10) ( (r/reff)^(1/n)-1 )
        #b_n taken based on solution to gamma function (Capaccioli 1989) for n < 10.
        return mueff + 2.5*(0.868*n-0.142)*((r/reff)**(1./n) - 1)

    def fit_sersic(self, r, surface_brightness, return_cov=False):
        s0_guess = np.mean(surface_brightness[:3])
        s0_range = [10,40]
        n_range = [0.5,16.5]
        r0_range=[0,100]

        r0_guess = min(r[int(len(r)/2)], 50)
        n_guess = 1.0

        sigma = 10**(0.6*(surface_brightness-20))/r
        sigma = None

        popt, pcov = scipy.optimize.curve_fit(self.sersic_surface_brightness,r,surface_brightness,
                                          bounds=np.array((s0_range, r0_range, n_range)).T,
                                          sigma=sigma,
                                          p0=(s0_guess, r0_guess, n_guess))

        if return_cov:
            return popt, pcov
        else:
            return popt

    def live_calculate(self, halo, *args):
        delta_r = self.get_simulation_property("approx_resolution_kpc", 0.1)
        r0 = delta_r/2
        sb_property = self.band+"_surface_brightness"
        if sb_property in halo.keys():
            surface_brightness = halo[sb_property]
        else:
            try:
                surface_brightness = halo.calculate(sb_property+"()")
            except NoResultsError:
                return null_result(self)
        if len(surface_brightness)<self.smooth*2:
            return null_result(self)
        if surface_brightness.max()==0:
            return null_result(self)
        r = np.arange(len(surface_brightness))*delta_r + delta_r/2.
        nbins = len(surface_brightness)
        if self.smooth > 1:
            surface_brightness_new = np.nanmean(np.pad(surface_brightness[self.smooth:].astype(float),
                                                (0,3-surface_brightness[self.smooth:].size%self.smooth),mode='constant',
                                                constant_values=np.NaN).reshape(-1,self.smooth),axis=1)
            surface_brightness = np.insert(surface_brightness_new,0,np.mean(surface_brightness[:self.smooth]))
            r = np.arange(len(surface_brightness))*(delta_r*self.smooth) + delta_r*self.smooth/2.
            r[0] = r0


        flux_density = 10**(surface_brightness/-2.5)
        flux_density[flux_density!=flux_density]=0
        cumu_flux_density = (r * flux_density).cumsum()
        cumu_flux_density/=cumu_flux_density[-1]

        try:
            half_light_i = np.where(cumu_flux_density>0.5)[0][0]
            half_light = r0+delta_r * half_light_i * self.smooth
        except:
            half_light = None
        if self.type=='hlr':
            if half_light is None:
                return null_result(self)
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
            return null_result(self)
        else:
            try:
                m0, r0, n = self.fit_sersic(r_fit, sb_fit)
            except:
                m0 = None
                n = None
                r0 = None
            return half_light, m0, n, r0


class VelDispersionProfileBase(BinnedProfileProperty):
    '''
    Shared machinery for the velocity dispersion profiles: subclasses differ
    only in how the variance of a velocity component is taken from a profile.
    '''
    names = None

    # Dispersions are noise below this many particles of a given species.
    min_particles = MIN_PARTICLES

    def plot_xlabel(self):
        return "r/kpc"

    def requires_property(self):
        return ["shrink_center", "max_radius"]

    def variance(self, pro, component):
        '''Mass weighted variance of the named velocity component in each bin'''
        raise NotImplementedError

    def dispersions(self, pro):
        '''Radial and total 3d velocity dispersion of one species' profile'''
        if pro is None:
            return None, None
        vx2, vy2, vz2 = (self.variance(pro, component) for component in ('vx', 'vy', 'vz'))
        return np.sqrt(self.variance(pro, 'vr')), np.sqrt(vx2+vy2+vz2)

    def rstat(self, halo, maxrad):
        delta, nbins, maxrad = self.binning(maxrad)
        #create new variables for the squares of the velocities, which the
        #profile then averages to give the mean square velocity in each bin
        for component in ('vr', 'vx', 'vy', 'vz'):
            halo[component+'2'] = halo[component]**2
        profiles = [family_profile(family, maxrad, nbins, min_particles=self.min_particles)
                    for family in (halo.s, halo.g, halo.dm)]
        dispersions = [self.dispersions(pro) for pro in profiles]
        return tuple(radial for radial, _ in dispersions) \
            + tuple(total for _, total in dispersions)

    @centred_calculation
    @velocity_centred()
    def calculate(self, halo, properties):
        return tuple(nan_remover(self.rstat(halo, properties['max_radius'])))


class VelDispersionProfile(VelDispersionProfileBase):
    '''
    dispersion of stars, gas and dm profiles
    '''
    names = "vrdisp_star", "vrdisp_gas", "vrdisp_dm", "vdisp_star_3d", "vdisp_gas_3d", "vdisp_dm_3d"

    def variance(self, pro, component):
        '''Mass weighted variance within each radial bin'''
        return pro[component+'2'] - pro[component]**2


class VelDispersionProfileEncl(VelDispersionProfileBase):
    '''
    enclosed velosity dispersion profiles
    '''
    names = "vrdisp_encl_star", "vrdisp_encl_gas", "vrdisp_encl_dm", "vdisp_encl_star_3d", "vdisp_encl_gas_3d", "vdisp_encl_dm_3d"

    min_particles = 5

    @staticmethod
    def enclosed_mean(pro, key):
        '''Weighted mean of a binned quantity over everything interior to each bin'''
        return np.nancumsum(pro[key]*pro['weight_fn'])/np.nancumsum(pro['weight_fn'])

    def variance(self, pro, component):
        '''Mass weighted variance of everything interior to each radial bin'''
        return self.enclosed_mean(pro, component+'2') \
            - self.enclosed_mean(pro, component)**2


#Define two profile properties, computing the enclosed dispersions directly from
#the particle data rather than from the binned means
def cumulative_variance(pro, key):
    variance = np.zeros(pro.nbins)
    cumind = []
    values = pro.sim[key].view(np.ndarray)
    weight = pro.sim[pro._weight_by].view(np.ndarray)
    for i in range(pro.nbins):
        cumind = np.append(cumind,pro.binind[i]).astype(np.int64)
        total_weight = weight[cumind].sum()
        mean_sq = ((values*weight)[cumind].sum()/total_weight)**2
        sq_mean = (values**2*weight)[cumind].sum()/total_weight
        variance[i] = sq_mean - mean_sq
    return variance

@pynbody.analysis.profile.Profile.profile_property
def vr_disp_encl(self):
    return np.sqrt(cumulative_variance(self, 'vr'))

@pynbody.analysis.profile.Profile.profile_property
def v_disp_tot_encl(self):
    return np.sqrt(sum(cumulative_variance(self, component)
                       for component in ('vx', 'vy', 'vz')))


class BHAccAveHistogram(BHAccHistogram):

    names = "BH_mdot_histogram_ave",

    def requires_property(self):
        return []


    def calculate(self, halo, properties):
        grid_tmax_Gyr = 20.0  # calculate up to 20 Gyr
        nbins = int(grid_tmax_Gyr / self.pixel_delta_t_Gyr)

        halo = halo.s

        if len(halo)!=1:
            raise RuntimeError("Not a BH!")

        if halo['tform'][0]>0:
            raise RuntimeError("Not a BH!")

        mask = self.log.vars['bhid']==halo['iord']
        if(mask.sum()==0):
            raise RuntimeError("Can't find BH in .orbit file")

        t_orbit = self.log.vars['time'][mask]
        Mdot_orbit = self.log.vars['mdotmean'][mask]
        order = np.argsort(t_orbit)

        t_max = properties.timestep.time_gyr

        mdot_grid_n, _ = np.histogram(t_orbit[order],bins=nbins,range=(0,grid_tmax_Gyr))
        mdot_grid_sum, _ = np.histogram(t_orbit[order],weights=Mdot_orbit[order],bins=nbins,range=(0,grid_tmax_Gyr))
        mdot_grid_ave = mdot_grid_sum/mdot_grid_n

        return mdot_grid_ave[self.store_slice(t_max)]
