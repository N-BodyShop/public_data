import os
import pytest
import tangos
import numpy as np
import pynbody as pyn
from numpy.testing import assert_allclose

REL_TOL=5e-2 # No calculations can produce results more than 5% different than what's in the database

@pytest.fixture(scope='function', params=os.listdir('testdata'))
def all_sims(request, pyn_snaps):
    tsim = tangos.get_simulation(request.param)
    snap = tsim.timesteps[0]
    db_h = tsim.timesteps[0].halos[0]
    sim = pyn_snaps[0][snap.filename]
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
    assert_allclose(cen, h['pos'].mean(axis=0).in_units('kpc'), rtol=REL_TOL)

def test_masses(all_sims):
    """
    Check that the tangos masses match the mass loaded straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    assert(db_h['finder_mass'] > 0)
    assert(db_h['finder_star_mass'] >= 0)
    assert(db_h['finder_gas_mass'] >= 0)
    assert(db_h['finder_dm_mass'] >= 0)
    assert(db_h['finder_mass'] == h['mass'].in_units('Msol').sum())
    assert(db_h['finder_star_mass'] == h.s['mass'].in_units('Msol').sum())
    assert(db_h['finder_gas_mass'] == h.g['mass'].in_units('Msol').sum())
    assert(db_h['finder_dm_mass'] == h.d['mass'].in_units('Msol').sum())

def test_spin_parameters(all_sims):
    """
    Check that the tangos spin parameters matches the values from AHF.
    """
    db_h, sim, h, cen, res = all_sims
    assert(db_h['lambda'] > 0)
    assert(db_h['lambda_gas'] > 0)
    assert(db_h['lambda_star'] > 0)
    assert(db_h['lambda'] == h.properties['lambda'].sum())
    assert(db_h['lambda_gas'] == h.properties['lambda_gas'].sum())
    assert(db_h['lambda_star'] == h.properties['lambda_star'].sum())

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
    if 'massHot' in sim.loadable_keys():
        eligible &= ~pyn.filt.HighPass('massHot', 0)
    elif 'MassHot' in sim.loadable_keys():
        eligible &= ~pyn.filt.HighPass('MassHot', 0)
    sfg = h.g[eligible]
    tdyn = 1.0/np.sqrt(4*np.pi*pyn.units.G*sfg['rho'])
    sfr = dCstar*(sfg['mass']/tdyn).in_units('Msol yr**-1')
    if H2:
        sfr *= sfg['H2']
    assert(sfr.sum() == db_h['instantaneous_SFR'])

def test_SFH(all_sims):
    """
    Check that the SFR matches the total stellar mass.
    """
    db_h, sim, h, cen, res = all_sims
    assert(np.all(db_h['SFR_histogram'] >= 0))
    # Only the last 500 Myr are stored by default with tangos when a single db_hshot is used.
    assert_allclose(db_h['SFR_histogram'].mean()*db_h.timestep.time_gyr*1e9, 
    h.s[h.s['tform'].in_units('Myr') > sim.properties['time'].in_units('Myr')-500]['massform'].in_units('Msol').sum(), rtol=REL_TOL)

def test_mass_profiles(all_sims):
    """
    Check that the mass profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    coldfilt = pyn.filt.LowPass('temp', 1.e5)
    warmfilt = pyn.filt.BandPass('temp', 1.e5, 1.e6)
    hotfilt = pyn.filt.HighPass('temp', 1.e6)
    if 'massHot' in h.loadable_keys():
        twophase = pyn.filt.HighPass('massHot', 0)
        sim.g['massCold'] = sim.g['mass'] - sim.g['massHot']
        hotmass = h.g[twophase]['massHot'].in_units('Msol').sum() 
    else:
        twophase = pyn.filt.HighPass('MassHot', 0)
        sim.g['massCold'] = sim.g['mass'] - sim.g['MassHot']
        hotmass = h.g[twophase]['MassHot'].in_units('Msol').sum() 
    hotmass += h.g[~twophase & hotfilt]['mass'].in_units('Msol').sum() 
    coldmass = h.g[~twophase & coldfilt]['mass'].in_units('Msol').sum() + h.g[twophase]['massCold'].in_units('Msol').sum() 
    warmmass = h.g[~twophase & warmfilt]['mass'].in_units('Msol').sum() 
    # Our gas temperature phases should be less than the total mass
    assert(np.all(db_h['gas_mass_profile'] >= db_h['cold_gas_mass_profile']))
    assert(np.all(db_h['gas_mass_profile'] >= db_h['warm_gas_mass_profile']))
    assert(np.all(db_h['gas_mass_profile'] >= db_h['hot_gas_mass_profile']))

    # The profile masses should match what's in the raw db_hshot
    assert_allclose(db_h['cold_gas_mass_profile'][-1],
    coldmass, rtol=REL_TOL)
    assert_allclose(db_h['warm_gas_mass_profile'][-1],
    warmmass, rtol=REL_TOL)
    assert_allclose(db_h['hot_gas_mass_profile'][-1],
    hotmass, rtol=REL_TOL)
    assert_allclose(db_h['gas_mass_profile'][-1],
    np.sum(h.g['mass'].in_units('Msol')), rtol=REL_TOL)
    assert_allclose(db_h['star_mass_profile'][-1],
    np.sum(h.s['mass'].in_units('Msol')), rtol=REL_TOL)
    assert_allclose(db_h['dm_mass_profile'][-1],
    np.sum(h.d['mass'].in_units('Msol')), rtol=REL_TOL)

def test_metal_profiles(all_sims):
    """
    Check that the metal fraction profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    twophase = pyn.filt.HighPass('massHot', 0)
    coldfilt = pyn.filt.LowPass('temp', 1.e5)
    for i in ['gas', 'cold_gas', 'star']:
        for j in ['metal', 'Fe', 'Ox']:
            # Fe/Ox Fraction and Metallicity is between zero and one
            assert(np.all(db_h[f'{i}_{j}_profile'] >= 0))
            assert(np.all(db_h[f'{i}_{j}_profile'] <= 1))
            # Metal fractions should always be greater than individual species fractions
            assert(np.all(db_h[f'{i}_metal_profile'] >= db_h[f'{i}_{j}_profile']))

    # Check that the metal fractions match the db_hshot totals.
    star_masses = db_h['star_mass_profile']
    star_masses[1:] -= star_masses[:-1]
    gas_masses = db_h['gas_mass_profile']
    gas_masses[1:] -= gas_masses[:-1]
    cold_gas_masses = db_h['cold_gas_mass_profile']
    cold_gas_masses[1:] -= cold_gas_masses[:-1]
    assert_allclose((h.s['mass'].in_units('Msol')*h.s['metals']).sum(),
    (star_masses*db_h['star_metal_profile']).sum(), rtol=REL_TOL)
    assert_allclose((h.g['mass'].in_units('Msol')*h.g['metals']).sum(),
    (gas_masses*db_h['gas_metal_profile']).sum(), rtol=REL_TOL)
    assert_allclose((h.s['mass'].in_units('Msol')*h.s['FeMassFrac']).sum(),
    (star_masses*db_h['star_Fe_profile']).sum(), rtol=REL_TOL)
    assert_allclose((h.g['mass'].in_units('Msol')*h.g['FeMassFrac']).sum(),
    (gas_masses*db_h['gas_Fe_profile']).sum(), rtol=REL_TOL)
    assert_allclose((h.s['mass'].in_units('Msol')*h.s['OxMassFrac']).sum(),
    (star_masses*db_h['star_Ox_profile']).sum(), rtol=REL_TOL)
    assert_allclose((h.g['mass'].in_units('Msol')*h.g['OxMassFrac']).sum(),
    (gas_masses*db_h['gas_Ox_profile']).sum(), rtol=REL_TOL)

def test_surface_brightness(all_sims):
    """
    Check that the stellar surface brightness profiles are sane
    """
    db_h, sim, h, cen, res = all_sims
    with pyn.analysis.angmom.faceon(h):
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
    assert_allclose((db_h['dm_density_profile']*shell_vol).sum(),
    h.d['mass'].in_units('Msol').sum(), rtol=REL_TOL)
    assert_allclose((db_h['star_density_profile']*shell_vol).sum(),
    h.s['mass'].in_units('Msol').sum(), rtol=REL_TOL)
    assert_allclose((db_h['gas_density_profile']*shell_vol).sum(),
    h.g['mass'].in_units('Msol').sum(), rtol=REL_TOL)

def test_inflow_outflow(all_sims):
    """
    Check that the inflow and outflow profiles match ones calculated straight from pynbody
    """
    db_h, sim, h, cen, res = all_sims
    with h.translate(-cen):
        with pyn.analysis.halo.vel_center(h, cen_size=db_h['max_radius']):
            # Check that all values are positive
            assert(np.all(db_h['mass_inflow_profile'] >= 0))
            assert(np.all(db_h['mass_outflow_profile'] >= 0))
            assert(np.all(db_h['metal_inflow_profile'] >= 0))
            assert(np.all(db_h['metal_outflow_profile'] >= 0))
            assert(np.all(db_h['Fe_inflow_profile'] >= 0))
            assert(np.all(db_h['Fe_outflow_profile'] >= 0))
            assert(np.all(db_h['Ox_inflow_profile'] >= 0))
            assert(np.all(db_h['Ox_outflow_profile'] >= 0))
            # Metal masses should never be larger than the total mass
            assert(np.all(db_h['mass_inflow_profile'] >= db_h['metal_inflow_profile']))
            assert(np.all(db_h['mass_outflow_profile'] >= db_h['metal_outflow_profile']))
            assert(np.all(db_h['metal_inflow_profile'] >= db_h['Fe_inflow_profile']))
            assert(np.all(db_h['metal_outflow_profile'] >= db_h['Fe_outflow_profile']))
            assert(np.all(db_h['metal_inflow_profile'] >= db_h['Ox_inflow_profile']))
            assert(np.all(db_h['metal_outflow_profile'] >= db_h['Ox_outflow_profile']))
            # Check that our totals match what is in the db_hshot
            iflow = h.g[h.g['vr'] < 0]
            oflow = h.g[h.g['vr'] > 0]
            assert_allclose(db_h['mass_inflow_profile'].sum(),
                                       (-iflow['vr'].in_units('kpc yr**-1')*iflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL)
            assert_allclose(db_h['mass_outflow_profile'].sum(),
                                       (oflow['vr'].in_units('kpc yr**-1')*oflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL)
            assert_allclose(db_h['metal_inflow_profile'].sum(),
                                       (-iflow['vr'].in_units('kpc yr**-1')*iflow['metals']*iflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 
            assert_allclose(db_h['metal_outflow_profile'].sum(),
                                       (oflow['vr'].in_units('kpc yr**-1')*oflow['metals']*oflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 
            assert_allclose(db_h['Fe_inflow_profile'].sum(),
                                       (-iflow['vr'].in_units('kpc yr**-1')*iflow['FeMassFrac']*iflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 
            assert_allclose(db_h['Fe_outflow_profile'].sum(),
                                       (oflow['vr'].in_units('kpc yr**-1')*oflow['FeMassFrac']*oflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 
            assert_allclose(db_h['Ox_inflow_profile'].sum(),
                                       (-iflow['vr'].in_units('kpc yr**-1')*iflow['OxMassFrac']*iflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 
            assert_allclose(db_h['Ox_outflow_profile'].sum(),
                                       (oflow['vr'].in_units('kpc yr**-1')*oflow['OxMassFrac']*oflow['mass'].in_units('Msol')).sum()/res, 
                                       rtol=REL_TOL) 

def test_vrdisp(all_sims):
    """
    Check that the velocity dispersion calculations are sane.
    """
    db_h, sim, h, cen, res = all_sims
    with h.translate(-cen):
        with pyn.analysis.halo.vel_center(h, cen_size='5 kpc'):
            filt = pyn.filt.LowPass('r', res)
            # The radial velocity dispersion must always be positive
            assert(np.all(db_h['vrdisp_stars'] >= 0))
            assert(np.all(db_h['vrdisp_gas'] >= 0))
            assert(np.all(db_h['vrdisp_dm'] >= 0))
            assert(np.all(db_h['vrdisp_encl_stars'] >= 0))
            assert(np.all(db_h['vrdisp_encl_gas'] >= 0))
            assert(np.all(db_h['vrdisp_encl_dm'] >= 0))
            assert(np.all(db_h['vdisp_stars_3d'] >= 0))
            assert(np.all(db_h['vdisp_gas_3d'] >= 0))
            assert(np.all(db_h['vdisp_dm_3d'] >= 0))
            assert(np.all(db_h['vdisp_encl_stars_3d'] >= 0))
            assert(np.all(db_h['vdisp_encl_gas_3d'] >= 0))
            assert(np.all(db_h['vdisp_encl_dm_3d'] >= 0))

            # Check that the dispersion results match the db_hshot
            assert_allclose(db_h['vrdisp_stars'][0],
            np.std(h[filt].s['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vrdisp_dm'][0],
            np.std(h[filt].d['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vrdisp_gas'][0],
            np.std(h[filt].g['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_stars_3d'][0],
            np.linalg.norm(np.array(np.std(h[filt].s['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_dm_3d'][0],
            np.linalg.norm(np.array(np.std(h[filt].d['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_gas_3d'][0],
            np.linalg.norm(np.array(np.std(h[filt].g['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
            # Check that the enclosed dispersion results match the db_hshot
            assert_allclose(db_h['vrdisp_encl_stars'][-1],
            np.std(h.s['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vrdisp_encl_dm'][-1],
            np.std(h.d['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vrdisp_encl_gas'][-1],
            np.std(h.g['vr'].in_units('km s**-1')), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_encl_stars_3d'][-1],
            np.linalg.norm(np.array(np.std(h.s['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_encl_dm_3d'][-1],
            np.linalg.norm(np.array(np.std(h.d['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
            assert_allclose(db_h['vdisp_encl_gas_3d'][-1],
            np.linalg.norm(np.array(np.std(h.g['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
