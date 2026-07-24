import os
import pytest
import tangos
import numpy as np
import pynbody as pyn
from numpy.testing import assert_allclose

REL_TOL=5e-2 # No calculations can produce results more than 5% different than what's in the database
auxnames = {'mass':'mass', 'metal':'metals',
            'Fe':'FeMassFrac','Ox':'OxMassFrac'} # Mappings to pynbody keys

def families(h):
    """
    return two zipped tuples of the three particle type names and their
    components from the provided halo.
    """
    return zip(['gas', 'star', 'dm'], [h.g, h.s, h.d])

@pytest.fixture(scope='function', params=['g3021','cptmarvel'])
def all_sims(request, pyn_snaps):
    tsim = tangos.get_simulation(request.param)
    snap = tsim.timesteps[0]
    db_h = tsim.timesteps[0].halos[0]
    sim = pyn_snaps[0][snap.filename]
    if 'MassHot' in sim.loadable_keys():
        sim.g['massHot'] = sim.g['MassHot'] 
    if 'massHot' not in sim.loadable_keys():
        sim.g['massHot'] = 0*sim.g['mass'] 
    sim.g['massCold'] = sim.g['mass'] - sim.g['massHot']
    h = pyn_snaps[1][snap.filename][0]
    cen = db_h['shrink_center']
    res = tsim['approx_resolution_kpc']
    return db_h, sim, h, cen, res

def test_max_radius(all_sims):
    """
    Check that the tangos max_radius matches the value loaded straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    with h.translate(-cen):
        assert_allclose(db_h['max_radius'], h['r'].max().in_units('kpc'), rtol=REL_TOL)

def test_center(all_sims):
    """
    Check that the tangos shrink_center matches the value loaded straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    assert_allclose(cen, pyn.analysis.halo.shrink_sphere_center(h.d, shrink_factor=0.8,
                                                                particles_for_velocity=0).in_units('kpc'), rtol=REL_TOL)

def test_masses(all_sims):
    """
    Check that the tangos masses match the mass loaded straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    assert(db_h[f'finder_mass'] >= 0)
    assert(db_h[f'finder_mass'] == h['mass'].in_units('Msol').sum())
    for i,j in families(h):
        assert(db_h[f'finder_{i}_mass'] >= 0)
        assert(db_h[f'finder_{i}_mass'] == j['mass'].in_units('Msol').sum())

def test_spin_parameters(all_sims):
    """
    Check that the tangos spin parameters matches the values from AHF.
    """
    db_h, sim, h, cen, res = all_sims
    for i in ['', '_gas', '_star']:
        assert(db_h[f'lambda{i}'] > 0)
        assert(db_h[f'lambda{i}'] == h.properties[f'lambda{i}'].sum())

def test_SFR(all_sims):
    """
    Check that the SFR matches the total stellar mass.
    """
    db_h, sim, h, cen, res = all_sims
    tsim = db_h.timestep.simulation
    dCstar = tsim.get('dCStar', 0.05)
    dPhysDenMin = tsim.get('dPhysDenMin',0.1)*pyn.units.m_p/(pyn.units.cm**3)
    dTempMax = tsim.get('dTempMax', 1.5e4)
    H2 = "COOLING_MOLECULARH" in tsim.get("macros", "")
    dense_filter = pyn.filt.HighPass('rho', dPhysDenMin)
    cold_filter = pyn.filt.LowPass('temp', dTempMax)
    eligible = dense_filter & cold_filter
    eligible &= ~pyn.filt.HighPass('massHot', 0)
    sfg = h.g[eligible]
    tdyn = 1.0/np.sqrt(4*np.pi*pyn.units.G*sfg['rho'])
    sfr = dCstar*(sfg['mass']/tdyn).in_units('Msol yr**-1')
    if H2:
        sfr *= sfg['H2']
    assert(sfr.sum() == db_h['instantaneous_SFR'])

def test_SFH(all_sims):
    """
    Check that the SFR histogram matches the star formation in the snapshot.
    """
    db_h, sim, h, cen, res = all_sims
    sfh = db_h['SFR_histogram']
    assert(np.all(sfh >= 0))
    # Tangos bins tform on a fixed 0-20 Gyr grid of histogram_delta_t_Gyr wide
    # bins, converts to Msol/yr, and stores only the final minimum_store_Gyr
    # (0.5 Gyr, i.e. the last nstored bins), zero padded back to t=0.  Rebuild
    # that same grid rather than cutting at t_now-500 Myr, which falls partway
    # through a bin and leaves the comparison good to only a few times 1e-4.
    dt = db_h.timestep.simulation.get('histogram_delta_t_Gyr')
    nstored = int(0.5/dt)
    start, end = len(sfh)-nstored, len(sfh)
    expected, _ = np.histogram(h.s['tform'].in_units('Gyr'),
                               weights=h.s['massform'].in_units('Msol'),
                               bins=int(20.0/dt), range=(0, 20.0))
    expected = expected/dt/1e9 # Msol per bin -> Msol/yr
    assert_allclose(sfh[start:end], expected[start:end], rtol=REL_TOL)
    # The stored bins must account for all the mass formed in that window
    tform = h.s['tform'].in_units('Gyr')
    window = (tform >= start*dt) & (tform < end*dt)
    assert_allclose(sfh.sum()*dt*1e9,
                    h.s['massform'].in_units('Msol')[window].sum(), rtol=REL_TOL)

def test_mass_profiles(all_sims):
    """
    Check that the mass profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    coldfilt = pyn.filt.LowPass('temp', 1.e5)
    warmfilt = pyn.filt.BandPass('temp', 1.e5, 1.e6)
    hotfilt = pyn.filt.HighPass('temp', 1.e6)
    mass = {'hot':0,'warm':0,'cold':0}
    twophase = pyn.filt.HighPass('massHot', 0)
    radfilt = pyn.filt.Sphere(db_h['max_radius'])
    with sim.translate(-cen):
        mass['hot'] = sim[radfilt].g[twophase]['massHot'].in_units('Msol').sum() 
        mass['hot'] += sim[radfilt].g[~twophase & hotfilt]['mass'].in_units('Msol').sum() 
        mass['cold'] = sim[radfilt].g[~twophase & coldfilt]['mass'].in_units('Msol').sum() \
        + h.g[twophase]['massCold'].in_units('Msol').sum() 
        mass['warm'] = sim[radfilt].g[~twophase & warmfilt]['mass'].in_units('Msol').sum() 
    for i in ['cold', 'warm', 'hot']:
        # Our gas temperature phases should be less than the total mass
        assert(np.all(db_h['gas_mass_profile'] >= db_h[f'{i}_gas_mass_profile']))

        # The profile masses should match what's in the raw snapshot
        assert_allclose(db_h[f'{i}_gas_mass_profile'][-1],
        mass[i], rtol=REL_TOL)
    for i,j in families(h):
        assert_allclose(db_h[f'{i}_mass_profile'][-1],
        np.sum(j['mass'].in_units('Msol')), rtol=REL_TOL)

def test_metal_profiles(all_sims):
    """
    Check that the metal fraction profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    twophase = pyn.filt.HighPass('massHot', 0)
    coldfilt = pyn.filt.LowPass('temp', 2.e4)
    with sim.translate(-cen):
        filt = pyn.filt.Sphere(db_h['max_radius'])
        for i,j in zip(['gas', 'cold_gas', 'star'], [sim[filt].g, sim[filt].g, sim[filt].s]):
            for k in ['metal', 'Fe', 'Ox']:
                # Fe/Ox Fraction and Metallicity is between zero and one
                assert(np.all(db_h[f'{i}_{k}_profile'] >= 0))
                assert(np.all(db_h[f'{i}_{k}_profile'] <= 1))
                # Metal fractions should always be greater than individual species fractions
                assert(np.all(db_h[f'{i}_metal_profile'] >= db_h[f'{i}_{k}_profile']))

                # Check that the metal fractions match the snapshot totals.
                if i != 'cold_gas':
                    masses = db_h[f'{i}_mass_profile']
                    masses[1:] -= masses[:-1]
                    assert_allclose((j['mass'].in_units('Msol')*j[auxnames[k]]).sum(),
                                    (masses*db_h[f'{i}_{k}_profile']).sum(),
                                    rtol=REL_TOL)

def test_surface_brightness(all_sims):
    """
    Check that the stellar surface brightness profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    with sim.translate(-cen):
        with pyn.analysis.halo.vel_center(h):
            with pyn.analysis.angmom.faceon(h,already_centered=True):
                rmax = db_h['max_radius']
                ps = pyn.analysis.profile.Profile(h.s, type='lin', ndim=2, rmin=0, rmax=rmax, nbins=int(rmax/res))
                pkeys = ['sb,'+x for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
                tkeys = [x+'_surface_brightness' for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
                for pk,tk in zip(pkeys,tkeys):
                    assert_allclose(db_h[tk], ps[pk], rtol=REL_TOL)

def test_density_profiles(all_sims):
    """
    Check that the mass density profiles match ones calculated straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    rmax = db_h['max_radius']
    shell_vol = np.power(np.arange(res,int(rmax/res)*res+res,res),3)*np.pi*4/3
    shell_vol[1:] -= shell_vol[:-1]
    for i,j in families(h):
        assert_allclose((db_h[f'{i}_density_profile']*shell_vol).sum(),
        j['mass'].in_units('Msol').sum(), rtol=REL_TOL)

def test_inflow_outflow(all_sims):
    """
    Check that the inflow and outflow profiles match ones calculated straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    with sim.translate(-cen):
        with pyn.analysis.halo.vel_center(sim, cen_size=db_h['max_radius']):
            filt = pyn.filt.Sphere(db_h['max_radius'])
            flow = {'inflow':sim[filt].g[sim[filt].g['vr'] < 0], 'outflow':sim[filt].g[sim[filt].g['vr'] > 0]}
            for i in ['mass', 'metal', 'Fe', 'Ox']:
                for j in ['inflow','outflow']:
                    # Check that all values are positive
                    assert(np.all(db_h[f'{i}_{j}_profile'] >= 0))
                    # Metal masses should never be larger than the total mass
                    assert(np.all(db_h[f'mass_{j}_profile'] >= db_h[f'{i}_{j}_profile']))
                    if i != 'mass':
                        assert(np.all(db_h[f'metal_{j}_profile'] >= db_h[f'{i}_{j}_profile']))
                    # Check that our totals match what is in the snapshot
                    weight = flow[j][auxnames[i]] if i != 'mass' else 1
                    assert_allclose(db_h[f'{i}_{j}_profile'].sum(), \
                            np.abs(flow[j]['vr'].in_units('kpc yr**-1') * 
                             weight*flow[j]['mass'].in_units('Msol')).sum()/res, rtol=REL_TOL)

def test_contamination(all_sims):
    """
    Check that the contamination fraction matches the heavy/light dark matter
    particle counts from the snapshot.
    """
    db_h, sim, h, cen, res = all_sims
    assert(db_h['contamination_fraction'] >= 0)
    assert(db_h['contamination_fraction'] <= 1)
    min_dm_mass = h.d['mass'].min().in_units('Msol')
    res_msol = db_h.timestep.simulation.get('approx_resolution_Msol', np.inf)
    if np.isclose(min_dm_mass, res_msol, rtol=1e-1, atol=1):
        expected = float((h.d['mass'].in_units('Msol') > min_dm_mass).sum())/len(h.d)
    else:
        # The halo contains no particles at the deepest zoom level
        expected = 1.0
    assert_allclose(db_h['contamination_fraction'], expected, rtol=REL_TOL)

def test_angmom_profiles(all_sims):
    """
    Check that the specific angular momentum profiles match ones calculated
    straight from pynbody. Profile properties are computed over all particles
    within a sphere of max_radius, not just halo members.
    """
    db_h, sim, h, cen, res = all_sims
    rmax = db_h['max_radius']
    with sim.translate(-cen):
        reg = sim[pyn.filt.Sphere(rmax)]
        with pyn.analysis.halo.vel_center(reg, cen_size='5 kpc'):
            for i,j in zip(['gas', 'star', 'dm'], [reg.g, reg.s, reg.d]):
                # Magnitudes are positive and angles fall in the expected
                # ranges (empty bins hold NaN, so compare only finite values)
                jmag = db_h[f'j_{i}_profile']
                jtheta = db_h[f'j_theta_{i}_profile']
                jphi = db_h[f'j_phi_{i}_profile']
                assert(np.all(jmag[np.isfinite(jmag)] >= 0))
                assert(np.all(jtheta[np.isfinite(jtheta)] >= 0))
                assert(np.all(jtheta[np.isfinite(jtheta)] <= np.pi))
                assert(np.all(np.abs(jphi[np.isfinite(jphi)]) <= np.pi))
                # Check that the profiles match the snapshot
                p = pyn.analysis.profile.Profile(j, type='lin', ndim=3,
                                                 min=0, max=rmax, nbins=int(rmax/res))
                assert_allclose(db_h[f'j_{i}_profile'], p['jtot'], rtol=REL_TOL)
                assert_allclose(db_h[f'j_theta_{i}_profile'], p['j_theta'], rtol=REL_TOL)
                assert_allclose(db_h[f'j_phi_{i}_profile'], p['j_phi'], rtol=REL_TOL)

def test_vrdisp(all_sims):
    """
    Check that the velocity dispersion calculations are sane.
    """
    db_h, sim, h, cen, res = all_sims
    normify = lambda x: np.linalg.norm(np.array(np.std(x, axis=0)))
    with sim.translate(-cen):
        with pyn.analysis.halo.vel_center(h, cen_size='5 kpc'):
            filt = pyn.filt.Sphere(res)
            # The radial velocity dispersion must always be positive
            for i,j in families(sim[pyn.filt.Sphere(db_h['max_radius'])]):
                for k in ['', 'encl_']:
                    assert(np.all(db_h[f'vrdisp_{k}{i}'] >= 0))
                    assert(np.all(db_h[f'vdisp_{k}{i}_3d'] >= 0))
                if len(j[filt]) == 0:
                    # Check that the dispersion results match the snapshot
                    assert_allclose(db_h[f'vrdisp_{i}'][0],
                                    0,
                                    rtol=REL_TOL)
                    assert_allclose(db_h[f'vdisp_{i}_3d'][0],
                                    0,
                                    rtol=REL_TOL)
                else:
                    # Check that the dispersion results match the snapshot
                    assert_allclose(db_h[f'vrdisp_{i}'][0],
                                    np.std(j[filt]['vr'].in_units('km s**-1')),
                                    rtol=REL_TOL)
                    assert_allclose(db_h[f'vdisp_{i}_3d'][0],
                                    normify(j[filt]['vel'].in_units('km s**-1')),
                                    rtol=REL_TOL)
                # Check that the enclosed dispersion results match the snapshot
                assert_allclose(db_h[f'vrdisp_encl_{i}'][-1],
                                np.std(j['vr'].in_units('km s**-1')),
                                rtol=REL_TOL)
                assert_allclose(db_h[f'vdisp_encl_{i}_3d'][-1],
                                normify(j['vel'].in_units('km s**-1')),
                                rtol=REL_TOL)
