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
        assert(snap.halos[0]['finder_mass'] == h[0]['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_star_mass'] == h[0].s['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_gas_mass'] == h[0].g['mass'].in_units('Msol').sum())
        assert(snap.halos[0]['finder_dm_mass'] == h[0].d['mass'].in_units('Msol').sum())
    assert(ran)