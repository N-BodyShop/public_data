# N-Body Shop Public Data Release Tools
[![Unit tests](https://github.com/N-BodyShop/public_data/actions/workflows/pytest.yml/badge.svg)](https://github.com/N-BodyShop/public_data/actions/workflows/pytest.yml)

The paper is started in [Overleaf](https://www.overleaf.com/project/66969b53942ad127acd619db) see @trquinn if you want to edit.

## Building the TANGOS Database
The simplest way to build a standard, N-Body Shop TANGOS database is to run the
`build_tangos_DB.sh CONFIGFILE` script.  The `CONFIGFILE` is a series of environment
variables that will specify what simulations to load, where the database should be
written, and any special options you want to use.  An example configuration file
is provided as `test.conf`.

### The Configuration File Format
The configuration file must contain the following variables:
- `TANGOS_SIMULATION_FOLDER`: The location containing your simulation outputs for TANGOS to read.
- `TANGOS_DB_CONNECTION`: The location where the TANGOS database should be written.
- `MIN_GAS`: The minimum number of gas particles a halo must have properties written in the database.
- `MIN_STARS`: The minimum number of star particles a halo must have properties written in the database.

The minimum number of gas/star particles can be set to zero if you want all halos to have properties written.
This can have a significant performance hit though, because many of the tangos/pynbody properties require
exceptions to be caught when there are no particles of a given type in a halo, and catching exceptions is expensive.
