import subprocess

def build(simname):
    return subprocess.run(['/bin/bash', '../build_tangos_DB.sh', simname, 'test.db'])

def test_build():
    build('/home/kellerbw/data/testdata').check_returncode()
