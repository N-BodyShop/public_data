import os
os.environ["TANGOS_SIMULATION_FOLDER"] = os.path.realpath("testdata")
os.environ["TANGOS_DB_CONNECTION"] = "test.db"
import tangos
import pynbody as pyn

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