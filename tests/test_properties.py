import os
import pytest
import tangos
import numpy as np
import pynbody as pyn

REL_TOL=5e-2 # No calculations can produce results more than 5% different than what's in the database

def test_max_radius(pyn_snaps):
    """
    Check that the tangos max_radius matches the value loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        cen = snap.halos[0]['shrink_center']
        with h[0].translate(-cen):
            np.testing.assert_allclose(snap.halos[0]['max_radius'], h[0]['r'].max().in_units('kpc'), rtol=REL_TOL)
    assert(ran)

def test_center(pyn_snaps):
    """
    Check that the tangos shrink_center matches the value loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        np.testing.assert_allclose(snap.halos[0]['shrink_center'],
        h[0]['pos'].mean(axis=0).in_units('kpc'), rtol=REL_TOL)
    assert(ran)

def test_masses(pyn_snaps):
    """
    Check that the tangos masses match the mass loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        assert(snap.halos[0]['finder_mass'] > 0)
        assert(snap.halos[0]['finder_star_mass'] >= 0)
        assert(snap.halos[0]['finder_gas_mass'] >= 0)
        assert(snap.halos[0]['finder_dm_mass'] >= 0)
        assert(snap.halos[0]['finder_mass'] == h[0]['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_star_mass'] == h[0].s['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_gas_mass'] == h[0].g['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_dm_mass'] == h[0].d['mass'].in_units('Msol').sum())
    assert(ran)

def test_spin_parameters(pyn_snaps):
    """
    Check that the tangos spin parameters matches the values from AHF.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        assert(snap.halos[0]['lambda'] > 0)
        assert(snap.halos[0]['lambda_gas'] > 0)
        assert(snap.halos[0]['lambda_star'] > 0)
        assert(snap.halos[0]['lambda'] == h[0].properties['lambda'].sum())
        assert(snap.halos[0]['lambda_gas'] == h[0].properties['lambda_gas'].sum())
        assert(snap.halos[0]['lambda_star'] == h[0].properties['lambda_star'].sum())
    assert(ran)

def test_SFR(pyn_snaps):
    """
    Check that the SFR matches the total stellar mass.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
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
        sfg = h[0].g[eligible]
        tdyn = 1.0/np.sqrt(4*np.pi*pyn.units.G*sfg['rho'])
        sfr = dCstar*(sfg['mass']/tdyn).in_units('Msol yr**-1')
        if H2:
            sfr *= sfg['H2']
        assert(sfr.sum() == snap.halos[0]['instantaneous_SFR'])
    assert(ran)

def test_SFH(pyn_snaps):
    """
    Check that the SFR matches the total stellar mass.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        assert(np.all(snap.halos[0]['SFR_histogram'] >= 0))
        # Only the last 500 Myr are stored by default with tangos when a single snapshot is used.
        np.testing.assert_allclose(snap.halos[0]['SFR_histogram'].mean()*snap.time_gyr*1e9, 
        h[0].s[h[0].s['tform'].in_units('Myr') > sim.properties['time'].in_units('Myr')-500]['massform'].in_units('Msol').sum(), rtol=REL_TOL)
    assert(ran)

def test_mass_profiles(pyn_snaps):
    """
    Check that the mass profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        coldfilt = pyn.filt.LowPass('temp', 1.e5)
        warmfilt = pyn.filt.BandPass('temp', 1.e5, 1.e6)
        hotfilt = pyn.filt.HighPass('temp', 1.e6)
        if 'massHot' in h[0].loadable_keys():
            twophase = pyn.filt.HighPass('massHot', 0)
            sim.g['massCold'] = sim.g['mass'] - sim.g['massHot']
            hotmass = h[0].g[twophase]['massHot'].in_units('Msol').sum() 
        else:
            twophase = pyn.filt.HighPass('MassHot', 0)
            sim.g['massCold'] = sim.g['mass'] - sim.g['MassHot']
            hotmass = h[0].g[twophase]['MassHot'].in_units('Msol').sum() 
        hotmass += h[0].g[~twophase & hotfilt]['mass'].in_units('Msol').sum() 
        coldmass = h[0].g[~twophase & coldfilt]['mass'].in_units('Msol').sum() + h[0].g[twophase]['massCold'].in_units('Msol').sum() 
        warmmass = h[0].g[~twophase & warmfilt]['mass'].in_units('Msol').sum() 
        # Our gas temperature phases should be less than the total mass
        assert(np.all(snap.halos[0]['gas_mass_profile'] >= snap.halos[0]['cold_gas_mass_profile']))
        assert(np.all(snap.halos[0]['gas_mass_profile'] >= snap.halos[0]['warm_gas_mass_profile']))
        assert(np.all(snap.halos[0]['gas_mass_profile'] >= snap.halos[0]['hot_gas_mass_profile']))

        # The profile masses should match what's in the raw snapshot
        np.testing.assert_allclose(snap.halos[0]['cold_gas_mass_profile'][-1],
        coldmass, rtol=REL_TOL)
        np.testing.assert_allclose(snap.halos[0]['warm_gas_mass_profile'][-1],
        warmmass, rtol=REL_TOL)
        np.testing.assert_allclose(snap.halos[0]['hot_gas_mass_profile'][-1],
        hotmass, rtol=REL_TOL)
        np.testing.assert_allclose(snap.halos[0]['gas_mass_profile'][-1],
        np.sum(h[0].g['mass'].in_units('Msol')), rtol=REL_TOL)
        np.testing.assert_allclose(snap.halos[0]['star_mass_profile'][-1],
        np.sum(h[0].s['mass'].in_units('Msol')), rtol=REL_TOL)
        np.testing.assert_allclose(snap.halos[0]['dm_mass_profile'][-1],
        np.sum(h[0].d['mass'].in_units('Msol')), rtol=REL_TOL)
    assert(ran)

def test_metal_profiles(pyn_snaps):
    """
    Check that the metal fraction profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        twophase = pyn.filt.HighPass('massHot', 0)
        coldfilt = pyn.filt.LowPass('temp', 1.e5)
        # Metallicity is between zero and one
        assert(np.all(snap.halos[0]['gas_metal_profile'] >= 0))
        assert(np.all(snap.halos[0]['gas_metal_profile'] <= 1))
        assert(np.all(snap.halos[0]['cold_gas_metal_profile'] >= 0))
        assert(np.all(snap.halos[0]['cold_gas_metal_profile'] <= 1))
        assert(np.all(snap.halos[0]['star_metal_profile'] >= 0))
        assert(np.all(snap.halos[0]['star_metal_profile'] <= 1))
        # Fe and Ox fraction is between zero and one
        assert(np.all(snap.halos[0]['gas_Fe_profile'] >= 0))
        assert(np.all(snap.halos[0]['gas_Fe_profile'] <= 1))
        assert(np.all(snap.halos[0]['cold_gas_Fe_profile'] >= 0))
        assert(np.all(snap.halos[0]['cold_gas_Fe_profile'] <= 1))
        assert(np.all(snap.halos[0]['star_Fe_profile'] >= 0))
        assert(np.all(snap.halos[0]['star_Fe_profile'] <= 1))
        assert(np.all(snap.halos[0]['gas_Ox_profile'] >= 0))
        assert(np.all(snap.halos[0]['gas_Ox_profile'] <= 1))
        assert(np.all(snap.halos[0]['cold_gas_Ox_profile'] >= 0))
        assert(np.all(snap.halos[0]['cold_gas_Ox_profile'] <= 1))
        assert(np.all(snap.halos[0]['star_Ox_profile'] >= 0))
        assert(np.all(snap.halos[0]['star_Ox_profile'] <= 1))
        # Metal fractions should always be greater than individual species fractions
        assert(np.all(snap.halos[0]['gas_metal_profile'] >= snap.halos[0]['gas_Fe_profile']))
        assert(np.all(snap.halos[0]['gas_metal_profile'] >= snap.halos[0]['gas_Ox_profile']))
        assert(np.all(snap.halos[0]['cold_gas_metal_profile'] >= snap.halos[0]['cold_gas_Fe_profile']))
        assert(np.all(snap.halos[0]['cold_gas_metal_profile'] >= snap.halos[0]['cold_gas_Ox_profile']))
        assert(np.all(snap.halos[0]['star_metal_profile'] >= snap.halos[0]['star_Fe_profile']))
        assert(np.all(snap.halos[0]['star_metal_profile'] >= snap.halos[0]['star_Ox_profile']))

        # Check that the metal fractions match the snapshot totals.
        star_masses = snap.halos[0]['star_mass_profile']
        star_masses[1:] -= star_masses[:-1]
        gas_masses = snap.halos[0]['gas_mass_profile']
        gas_masses[1:] -= gas_masses[:-1]
        cold_gas_masses = snap.halos[0]['cold_gas_mass_profile']
        cold_gas_masses[1:] -= cold_gas_masses[:-1]
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['metals']).sum(),
        (star_masses*snap.halos[0]['star_metal_profile']).sum(), rtol=REL_TOL)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['metals']).sum(),
        (gas_masses*snap.halos[0]['gas_metal_profile']).sum(), rtol=REL_TOL)
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['FeMassFrac']).sum(),
        (star_masses*snap.halos[0]['star_Fe_profile']).sum(), rtol=REL_TOL)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['FeMassFrac']).sum(),
        (gas_masses*snap.halos[0]['gas_Fe_profile']).sum(), rtol=REL_TOL)
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['OxMassFrac']).sum(),
        (star_masses*snap.halos[0]['star_Ox_profile']).sum(), rtol=REL_TOL)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['OxMassFrac']).sum(),
        (gas_masses*snap.halos[0]['gas_Ox_profile']).sum(), rtol=REL_TOL)
    assert(ran)

def test_surface_brightness(pyn_snaps):
    """
    Check that the stellar surface brightness profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        with pyn.analysis.angmom.faceon(h[0]):
            res = tsim['approx_resolution_kpc']
            rmax = snap.halos[0]['max_radius']
            ps = pyn.analysis.profile.Profile(h[0].s, type='lin', ndim=2, rmin=0, rmax=rmax, nbins=int(rmax/res))
            pkeys = ['sb,'+x for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
            tkeys = [x+'_surface_brightness' for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
            for pk,tk in zip(pkeys,tkeys):
                np.testing.assert_allclose(snap.halos[0][tk], ps[pk], rtol=REL_TOL)
    assert(ran)

def test_density_profiles(pyn_snaps):
    """
    Check that the mass density profiles match ones calculated straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        res = tsim['approx_resolution_kpc']
        rmax = snap.halos[0]['max_radius']
        shell_vol = np.power(np.arange(res,int(rmax/res)*res+res,res),3)*np.pi*4/3
        shell_vol[1:] -= shell_vol[:-1]
        np.testing.assert_allclose((snap.halos[0]['dm_density_profile']*shell_vol).sum(),
        h[0].d['mass'].in_units('Msol').sum(), rtol=REL_TOL)
        np.testing.assert_allclose((snap.halos[0]['star_density_profile']*shell_vol).sum(),
        h[0].s['mass'].in_units('Msol').sum(), rtol=REL_TOL)
        np.testing.assert_allclose((snap.halos[0]['gas_density_profile']*shell_vol).sum(),
        h[0].g['mass'].in_units('Msol').sum(), rtol=REL_TOL)
    assert(ran)

def test_inflow_outflow(pyn_snaps):
    """
    Check that the inflow and outflow profiles match ones calculated straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        cen = snap.halos[0]['shrink_center']
        with h[0].translate(-cen):
            with pyn.analysis.halo.vel_center(h[0], cen_size=snap.halos[0]['max_radius']):
                res = tsim['approx_resolution_kpc']
                # Check that all values are positive
                assert(np.all(snap.halos[0]['mass_inflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['mass_outflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['metal_inflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['metal_outflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['Fe_inflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['Fe_outflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['Ox_inflow_profile'] >= 0))
                assert(np.all(snap.halos[0]['Ox_outflow_profile'] >= 0))
                # Metal masses should never be larger than the total mass
                assert(np.all(snap.halos[0]['mass_inflow_profile'] >= snap.halos[0]['metal_inflow_profile']))
                assert(np.all(snap.halos[0]['mass_outflow_profile'] >= snap.halos[0]['metal_outflow_profile']))
                assert(np.all(snap.halos[0]['metal_inflow_profile'] >= snap.halos[0]['Fe_inflow_profile']))
                assert(np.all(snap.halos[0]['metal_outflow_profile'] >= snap.halos[0]['Fe_outflow_profile']))
                assert(np.all(snap.halos[0]['metal_inflow_profile'] >= snap.halos[0]['Ox_inflow_profile']))
                assert(np.all(snap.halos[0]['metal_outflow_profile'] >= snap.halos[0]['Ox_outflow_profile']))
                # Check that our totals match what is in the snapshot
                iflow = h[0].g[h[0].g['vr'] < 0]
                oflow = h[0].g[h[0].g['vr'] > 0]
                np.testing.assert_allclose(snap.halos[0]['mass_inflow_profile'].sum(),
                                           (-iflow['vr'].in_units('kpc yr**-1')*iflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['mass_outflow_profile'].sum(),
                                           (oflow['vr'].in_units('kpc yr**-1')*oflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['metal_inflow_profile'].sum(),
                                           (-iflow['vr'].in_units('kpc yr**-1')*iflow['metals']*iflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
                np.testing.assert_allclose(snap.halos[0]['metal_outflow_profile'].sum(),
                                           (oflow['vr'].in_units('kpc yr**-1')*oflow['metals']*oflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
                np.testing.assert_allclose(snap.halos[0]['Fe_inflow_profile'].sum(),
                                           (-iflow['vr'].in_units('kpc yr**-1')*iflow['FeMassFrac']*iflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
                np.testing.assert_allclose(snap.halos[0]['Fe_outflow_profile'].sum(),
                                           (oflow['vr'].in_units('kpc yr**-1')*oflow['FeMassFrac']*oflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
                np.testing.assert_allclose(snap.halos[0]['Ox_inflow_profile'].sum(),
                                           (-iflow['vr'].in_units('kpc yr**-1')*iflow['OxMassFrac']*iflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
                np.testing.assert_allclose(snap.halos[0]['Ox_outflow_profile'].sum(),
                                           (oflow['vr'].in_units('kpc yr**-1')*oflow['OxMassFrac']*oflow['mass'].in_units('Msol')).sum()/res, 
                                           rtol=REL_TOL) 
    assert(ran)

def test_vrdisp(pyn_snaps):
    """
    Check that the velocity dispersion calculations are sane.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn_snaps[0][snap.filename]
        h = pyn_snaps[1][snap.filename]
        res = tsim['approx_resolution_kpc']
        cen = snap.halos[0]['shrink_center']
        with h[0].translate(-cen):
            with pyn.analysis.halo.vel_center(h[0], cen_size='5 kpc'):
                filt = pyn.filt.LowPass('r', res)
                # The radial velocity dispersion must always be positive
                assert(np.all(snap.halos[0]['vrdisp_stars'] >= 0))
                assert(np.all(snap.halos[0]['vrdisp_gas'] >= 0))
                assert(np.all(snap.halos[0]['vrdisp_dm'] >= 0))
                assert(np.all(snap.halos[0]['vrdisp_encl_stars'] >= 0))
                assert(np.all(snap.halos[0]['vrdisp_encl_gas'] >= 0))
                assert(np.all(snap.halos[0]['vrdisp_encl_dm'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_stars_3d'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_gas_3d'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_dm_3d'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_encl_stars_3d'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_encl_gas_3d'] >= 0))
                assert(np.all(snap.halos[0]['vdisp_encl_dm_3d'] >= 0))

                # Check that the dispersion results match the snapshot
                np.testing.assert_allclose(snap.halos[0]['vrdisp_stars'][0],
                np.std(h[0][filt].s['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vrdisp_dm'][0],
                np.std(h[0][filt].d['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vrdisp_gas'][0],
                np.std(h[0][filt].g['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_stars_3d'][0],
                np.linalg.norm(np.array(np.std(h[0][filt].s['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_dm_3d'][0],
                np.linalg.norm(np.array(np.std(h[0][filt].d['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_gas_3d'][0],
                np.linalg.norm(np.array(np.std(h[0][filt].g['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
                # Check that the enclosed dispersion results match the snapshot
                np.testing.assert_allclose(snap.halos[0]['vrdisp_encl_stars'][-1],
                np.std(h[0].s['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vrdisp_encl_dm'][-1],
                np.std(h[0].d['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vrdisp_encl_gas'][-1],
                np.std(h[0].g['vr'].in_units('km s**-1')), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_encl_stars_3d'][-1],
                np.linalg.norm(np.array(np.std(h[0].s['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_encl_dm_3d'][-1],
                np.linalg.norm(np.array(np.std(h[0].d['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)
                np.testing.assert_allclose(snap.halos[0]['vdisp_encl_gas_3d'][-1],
                np.linalg.norm(np.array(np.std(h[0].g['vel'].in_units('km s**-1'),axis=0))), rtol=REL_TOL)

    assert(ran)
