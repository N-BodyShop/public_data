import tarfile
import requests
import subprocess
import pytest

def build(conf):
    return subprocess.run(['/bin/bash', 'build_tangos_DB.sh', conf])

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

@pytest.mark.order("first")
def test_build():
    get_testdata()
    build('test.conf').check_returncode()
