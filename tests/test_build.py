import pytest

def test_build(build_database):
    retcode = build_database
    retcode.check_returncode()
