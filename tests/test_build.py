import os
import pytest
import subprocess

def build(simname):
    return subprocess.run(['bash', '../build_tangos_DB.sh', simname, 'test.db'])

def test_build():
    os.environ["PYTHONPATH"] = '..'
    assert(build('testdata') == 0)
