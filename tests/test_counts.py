import os
os.environ["TANGOS_DB_CONNECTION"] = "test.db"
import tangos

def test_sim_count():
    """
    Make sure that the correct number of simulations are loaded.
    """
    assert(len(tangos.all_simulations()) == 2)

def test_snap_count():
    """
    Make sure that the correct number of snapshots are loaded.
    """
    assert(len(tangos.get_simulation('g3021').timesteps) == 1)
    assert(len(tangos.get_simulation('cptmarvel').timesteps) == 1)

def test_halo_count():
    """
    Make sure we have the correct number of halos.
    """
    assert(tangos.get_simulation('g3021').timesteps[0].halos.count() == 68)
    assert(tangos.get_simulation('cptmarvel').timesteps[0].halos.count() == 949)

def test_property_count():
    """
    Make sure we have the correct number of properties
    """
    N_import = len([line for line in open('import_properties').readlines() if line.strip() and line[0] != '#'])
    N_write = len([line for line in open('write_properties').readlines() if line.strip() and line[0] != '#'])
    assert(len(tangos.get_simulation('g3021').timesteps[0].halos[0].keys()) == N_import + N_write)
    assert(len(tangos.get_simulation('cptmarvel').timesteps[0].halos[0].keys()) == N_import + N_write)
