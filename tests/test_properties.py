import os
os.environ["TANGOS_SIMULATION_FOLDER"] = os.path.realpath("testdata")
os.environ["TANGOS_DB_CONNECTION"] = "test.db"
import tangos
import numpy as np
import pynbody as pyn

def test_max_radius():
    """
    Check that the tangos max_radius matches the value loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        pyn.analysis.center(h[0])
        np.testing.assert_allclose(snap.halos[0]['max_radius'], h[0]['r'].max().in_units('kpc'), rtol=1e-2)
    assert(ran)

def test_center():
    """
    Check that the tangos shrink_center matches the value loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        np.testing.assert_allclose(snap.halos[0]['shrink_center'],
        h[0]['pos'].mean(axis=0).in_units('kpc'), rtol=1e-2)
    assert(ran)

def test_masses():
    """
    Check that the tangos masses match the mass loaded straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        assert(snap.halos[0]['finder_mass'] > 0)
        assert(snap.halos[0]['finder_star_mass'] >= 0)
        assert(snap.halos[0]['finder_gas_mass'] >= 0)
        assert(snap.halos[0]['finder_dm_mass'] >= 0)
        assert(snap.halos[0]['finder_mass'] == h[0]['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_star_mass'] == h[0].s['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_gas_mass'] == h[0].g['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_dm_mass'] == h[0].d['mass'].in_units('Msol').sum())
    assert(ran)

def test_spin_parameters():
    """
    Check that the tangos spin parameters matches the values from AHF.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        assert(snap.halos[0]['lambda'] > 0)
        assert(snap.halos[0]['lambda_gas'] > 0)
        assert(snap.halos[0]['lambda_star'] > 0)
        assert(snap.halos[0]['lambda'] == h[0].properties['lambda'].sum())
        assert(snap.halos[0]['lambda_gas'] == h[0].properties['lambda_gas'].sum())
        assert(snap.halos[0]['lambda_star'] == h[0].properties['lambda_star'].sum())
    assert(ran)

def test_SFR():
    """
    Check that the SFR matches the total stellar mass.
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        assert(np.all(snap.halos[0]['SFR_histogram'] >= 0))
        SFR_hist_mass = snap.halos[0]['SFR_histogram'].mean()*snap.time_gyr*1e9
        snap_massform = h[0].s['massform'].in_units('Msol').sum()
        assert(SFR_hist_mass == snap_massform)
    assert(ran)

def test_mass_profiles():
    """
    Check that the mass profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
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
        coldmass, rtol=1e-2)
        np.testing.assert_allclose(snap.halos[0]['warm_gas_mass_profile'][-1],
        warmmass, rtol=1e-2)
        np.testing.assert_allclose(snap.halos[0]['hot_gas_mass_profile'][-1],
        hotmass, rtol=1e-2)
        np.testing.assert_allclose(snap.halos[0]['gas_mass_profile'][-1],
        np.sum(h[0].g['mass'].in_units('Msol')), rtol=1e-2)
        np.testing.assert_allclose(snap.halos[0]['star_mass_profile'][-1],
        np.sum(h[0].s['mass'].in_units('Msol')), rtol=1e-2)
        np.testing.assert_allclose(snap.halos[0]['dm_mass_profile'][-1],
        np.sum(h[0].d['mass'].in_units('Msol')), rtol=1e-2)
    assert(ran)

def test_metal_profiles():
    """
    Check that the metal fraction profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        twophase = pyn.filt.HighPass('massHot', 0)
        coldfilt = pyn.filt.LowPass('temp', 1.e5)
        sim.g['massCold'] = sim.g['mass'] - sim.g['massHot']
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

        # Check that the metal fractions match the snapshot totals.
        star_masses = snap.halos[0]['star_mass_profile']
        star_masses[1:] -= star_masses[:-1]
        gas_masses = snap.halos[0]['gas_mass_profile']
        gas_masses[1:] -= gas_masses[:-1]
        cold_gas_masses = snap.halos[0]['cold_gas_mass_profile']
        cold_gas_masses[1:] -= cold_gas_masses[:-1]
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['metals']).sum(),
        (star_masses*snap.halos[0]['star_metal_profile']).sum(), rtol=1e-2)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['metals']).sum(),
        (gas_masses*snap.halos[0]['gas_metal_profile']).sum(), rtol=1e-2)
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['FeMassFrac']).sum(),
        (star_masses*snap.halos[0]['star_Fe_profile']).sum(), rtol=1e-2)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['FeMassFrac']).sum(),
        (gas_masses*snap.halos[0]['gas_Fe_profile']).sum(), rtol=1e-2)
        np.testing.assert_allclose((h[0].s['mass'].in_units('Msol')*h[0].s['OxMassFrac']).sum(),
        (star_masses*snap.halos[0]['star_Ox_profile']).sum(), rtol=1e-2)
        np.testing.assert_allclose((h[0].g['mass'].in_units('Msol')*h[0].g['OxMassFrac']).sum(),
        (gas_masses*snap.halos[0]['gas_Ox_profile']).sum(), rtol=1e-2)
    assert(ran)

def test_surface_brightness():
    """
    Check that the stellar surface brightness profiles are sane
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        pyn.analysis.angmom.faceon(h[0])
        sim.physical_units()
        res = tsim['approx_resolution_kpc']
        rmax = snap.halos[0]['max_radius']
        ps = pyn.analysis.profile.Profile(h[0].s, type='lin', ndim=2, min=0, max=rmax, nbins=int(rmax/res))
        pkeys = ['sb,'+x for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
        tkeys = [x+'_surface_brightness' for x in ('u','g','r','i', 'z', 'U', 'V', 'J')]
        for pk,tk in zip(pkeys,tkeys):
            np.testing.assert_allclose(snap.halos[0][tk], ps[pk], rtol=5e-2)
    assert(ran)

def test_density_profiles():
    """
    Check that the mass density profiles match ones calculated straight from pynbody
    """
    ran = False
    for tsim in tangos.all_simulations():
        ran = True
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        h = sim.halos()
        res = tsim['approx_resolution_kpc']
        rmax = snap.halos[0]['max_radius']
        shell_vol = np.power(np.arange(res,int(rmax/res)*res+res,res),3)*np.pi*4/3
        shell_vol[1:] -= shell_vol[:-1]
        np.testing.assert_allclose((snap.halos[0]['dm_density_profile']*shell_vol).sum(),
        h[0].d['mass'].in_units('Msol').sum(), rtol=1e-2)
        np.testing.assert_allclose((snap.halos[0]['star_density_profile']*shell_vol).sum(),
        h[0].s['mass'].in_units('Msol').sum(), rtol=1e-2)
        np.testing.assert_allclose((snap.halos[0]['gas_density_profile']*shell_vol).sum(),
        h[0].g['mass'].in_units('Msol').sum(), rtol=1e-2)
    assert(ran)