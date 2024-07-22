#/usr/bin/env bash

# Import the data
#for i in $TANGOS_SIMULATION_FOLDER/MUGS2*
#do
#    tangos add `basename $i`
#done
tangos add MUGS2_g1536

# Set any non-default simulation properties for all imported simulations
./set_simulation_parameters.py

# Link to build the merger tree
##tangos link

# Import AHF Halo catalogue properties
##tangos import-properties `cat import_properties`

# Write galaxy properties
tangos write `cat write_properties`
