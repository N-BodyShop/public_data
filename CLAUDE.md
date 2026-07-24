# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Tools for building N-Body Shop public data release TANGOS databases from tipsy simulation snapshots (with AHF halo catalogues), using pynbody for the underlying analysis. There is no build step or package install; the repo is a set of scripts plus a tangos property-calculation module.

## Commands

```bash
# Build a TANGOS database (the main entry point)
./build_tangos_DB.sh CONFIGFILE        # e.g. test.conf, test_parallel.conf, MUGS2.conf

# Run tests (serial and parallel configs, as CI does), from the repo root.
# ALL FOUR variables must be set before pytest starts (see test notes below).
export TANGOS_DB_CONNECTION=test.db
export TANGOS_SIMULATION_FOLDER=$(realpath testdata)
export TANGOS_PROPERTY_MODULES=properties
export PYTHONPATH=.
pytest --paramfile=test.conf
pytest --paramfile=test_parallel.conf

# Run a single test
pytest tests/test_properties.py::test_masses --paramfile=test.conf
```

`--paramfile` is a custom option defined in `tests/conftest.py` (default `test.conf`); it selects the config the session-scoped build fixture passes to `build_tangos_DB.sh`. A `.venv` symlink to the working virtualenv exists in the repo root (gitignored).

Dependencies (installed via pip, no requirements file): `pynbody tangos scipy requests pytest pytest-order`.

Test notes:
- The session-scoped autouse fixtures in `tests/conftest.py` download ~test data from https://nbody.shop/testdata.tar.gz into `testdata/` and run `build_tangos_DB.sh` before any test runs, so even a single test triggers a full database build and a run takes many minutes end to end. `test.db` is deleted at session end.
- Tests expect to run from the repo root and need `TANGOS_SIMULATION_FOLDER` (absolute path to `testdata`), `TANGOS_DB_CONNECTION=test.db`, `TANGOS_PROPERTY_MODULES=properties`, and `PYTHONPATH=.` set in the shell **before pytest starts** (as CI does). The `os.environ` assignments in conftest/test modules are NOT sufficient: tangos reads `TANGOS_DB_CONNECTION` at import time during collection, so if it is unset when pytest launches, the test process silently queries the default `~/tangos_data.db` (empty) while the build subprocess correctly writes `test.db` — every test then fails with "No simulation matches" even though the build succeeded.
- If a previous pytest session was interrupted, a stale `test.db` may be left behind (the fixture that deletes it only runs at normal session end). Delete it before re-running so the build starts from a clean database rather than layering onto a stale, possibly partially-built one.
- `tests/test_counts.py` hardcodes the contents of the downloaded test data: 2 simulations (`g3021`, `cptmarvel`), 1 timestep each, 68 and 949 halos. If the tarball at nbody.shop changes, these fail for reasons unrelated to the code.

## Architecture

The database build pipeline (`build_tangos_DB.sh`) sources the config file as environment variables, then runs, in order:
1. `tangos add` for each subdirectory of `$TANGOS_SIMULATION_FOLDER`
2. `set_simulation_parameters.py` — stores per-simulation metadata (paramfile keys/values, raw logfile, log steps) on each simulation. These use a `_noweb` suffix so the tangos web view skips them.
3. `tangos link` — builds the merger tree
4. `tangos import-properties` — imports halo-finder (AHF) catalogue properties listed in the `import_properties` file
5. `tangos write` — computes properties listed in the `write_properties` file, filtered by `--include-only NGas()>$MIN_GAS` and `NStar()>$MIN_STAR`

Parallelism (`NPROCS`, `MPI`, `SERVER` config variables) maps to tangos backend flags: multiprocessing vs mpi4py, and `--load-mode=server`/`server-shared-mem` for volumes too large to load one snapshot per process.

Config-file gotchas:
- The script sources the config with `export $(cat $1 | grep -v '^#' | xargs)`, so values must not contain spaces or quotes; only `#` at line start comments out a line.
- The particle-count floors read by the script are `MIN_GAS` and `MIN_STAR`. The README documents the second as `MIN_STARS`, which the script ignores (it then defaults to 0).
- The script `unset PYTEST_CURRENT_TEST`s before running tangos: tangos refuses to load external property modules when that variable is set, so without it every custom property silently goes missing when the build runs under pytest.

`properties.py` is the custom tangos property module, activated via `TANGOS_PROPERTY_MODULES=properties` with the repo root on `PYTHONPATH` (the build script sets both). It defines:
- Profile-based properties subclassing `BinnedProfileProperty` (a `HaloDensityProfile` whose `binning()` gives every profile the same bins): inflow/outflow fluxes, angular momentum, metallicity, temperature-binned enclosed mass, surface brightness, velocity dispersions. The two dispersion classes share `VelDispersionProfileBase` and differ only in their `variance()`.
- Module-level helpers those classes are built from: `lin_profile`/`family_profile`/`profile_values` (profile construction, particle-count guards and result ordering), `velocity_centred` (decorator that removes and restores the halo bulk velocity), `two_phase_keys`/`hot_phase_key`, `load_metal_arrays`, `null_result` and `nan_remover`. Plus custom pynbody `Profile.profile_property` helpers: `p_r_*` (used by `InflowOutflow`) and `vr_disp_encl`/`v_disp_tot_encl` (registered but not currently used by any property).
- Live-calculated properties (`LiveRadius`, `StellarProfileDiagnosis` for sersic fits) that derive results from already-stored profiles at query time and are not listed in `write_properties`.
- `InstantaneousSFR`, which reads star-formation parameters (`dCStar`, `dPhysDenMin`, `dTempMax`, `macros`) from simulation properties via `get_simulation_property`.

Not every name in `write_properties` lives in `properties.py` — `contamination_fraction`, `shrink_center`, `max_radius`, `SFR_histogram`, the `*_density_profile`/`*_mass_profile` families and the `finder_*` masses are tangos built-ins. Conversely `properties.py` defines classes that are deliberately not written (`BHAccAveHistogram`, and the live-calculated ones).

To add a new stored property: implement it in `properties.py` (its `names` entries are the property names) and add those names to `write_properties`. Shared base classes must set `names = None`, otherwise tangos registers them as providers of whatever `names` they inherit (this is why `BinnedProfileProperty` and `VelDispersionProfileBase` do). Halo-finder properties go in `import_properties` instead. `tests/test_counts.py::test_property_count` asserts that halos have exactly the number of properties listed in those two files, so it will fail if a property is listed but fails to compute.

Domain conventions worth knowing:
- Two-phase gas snapshots store hot-phase mass under either `massHot` or `MassHot` (loader-dependent); `two_phase_keys` resolves the spelling and derives `massCold = mass - massHot`. Non-two-phase snapshots have neither, and code falls back to plain temperature cuts.
- Profile calculations return `None` per-array (not exceptions) when a halo lacks enough particles of a species; the particle-count guard in `family_profile` exists because catching pynbody exceptions for empty families is expensive.
- Velocity-centering calculations must restore the halo velocity after profiling, since the snapshot is shared across property calculations; `velocity_centred` does this in a `finally` block. A calculation that raises inside `@centred_calculation` also leaks the position re-centring (tangos' `_recenter` has no `finally`), corrupting later properties for that halo, so a failed `vel_center` returns nulls rather than propagating.
