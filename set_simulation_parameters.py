#!/usr/bin/env python3.12

import pynbody 
import tangos 
import numpy as np

# This short script is used to set all simulation
# parameters required by TANGOS to match the public
# data release choices.
if __name__ == "__main__":
    for sim in tangos.all_simulations():
        sim['histogram_delta_t_Gyr'] = 0.01
        sim['tangos_version'] = tangos.__version__
        sim['pynbody_version'] = pynbody.__version__
        # Load the parameter and log file
        ppath = sim.get_output_handler()._get_paramfile_path()
        pdict = sim.get_output_handler()._param_file_to_dict(ppath)
        lpath, prop_dict = sim.get_output_handler()._get_log_path(ppath, pdict)
        logfile = open(lpath)
        # The default Tangos ORM doesn't support dictionaries here, so we will 
        # just store keys and values as two lists
        sim['paramfile_keys'] = list(pdict.keys())
        sim['paramfile_vals'] = [pdict[key] for key in sim['paramfile_keys']]
        # Store the raw logfile as one long string
        sim['logfile'] = logfile.read()
        logfile.seek(0)
        # Store individual log steps as a numpy array.
        sim['logsteps'] = np.array([np.fromstring(line, sep=' ') for line in logfile.readlines() if 
                ('#' not in line and ('-' in line or '+' in line))])
