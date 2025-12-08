import os
import tangos
import pytest
import tarfile
import requests
import subprocess
import pynbody as pyn


@pytest.fixture(autouse=True)
def environment():
    os.environ["TANGOS_SIMULATION_FOLDER"] = os.path.realpath("testdata")
    os.environ["TANGOS_DB_CONNECTION"] = "test.db"

@pytest.fixture(scope="session")
def pyn_snaps():
    sims = {}
    halos = {}
    for tsim in tangos.all_simulations():
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        sim.physical_units()
        h = sim.halos()
        sims[snap.filename] = sim
        halos[snap.filename] = h
    return sims, halos

@pytest.fixture(scope="session")
def get_testdata():
    url = "https://nbody.shop/testdata.tar.gz"
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open("testdata.tar.gz", "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        # Extract the tar.gz file
        with tarfile.open("testdata.tar.gz", "r:gz") as tar:
            tar.extractall()
    else:
        raise Exception(f"Failed to download test data: {r.status_code}")

@pytest.fixture(scope="session", params=["test.conf", ], autouse=True)
def build_database(request, get_testdata):
    print(f"Building DB for {request.param}")
    yield subprocess.run(['/bin/bash', 'build_tangos_DB.sh', request.param])
    os.remove('test.db')
