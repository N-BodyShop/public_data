#!/usr/bin/env python3

import pynbody 
import tangos 

# This short script is used to set all simulation
# parameters required by TANGOS to match the public
# data release choices.
if __name__ == "__main__":
    for sim in tangos.all_simulations():
        sim['histogram_delta_t_Gyr'] = 0.01
        sim['tangos_version'] = tangos.__version__
        sim['pynbody_version'] = pynbody.__version__
