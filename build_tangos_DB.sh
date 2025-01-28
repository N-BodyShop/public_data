#/usr/bin/env bash

export TANGOS_SIMULATION_FOLDER=$1
export TANGOS_DB_CONNECTION=$2
export TANGOS_PROPERTY_MODULES=properties
export PYTHONPATH=.

# Import the data
for i in $TANGOS_SIMULATION_FOLDER/*
do
    tangos add `basename $i`
done

# Set any non-default simulation properties for all imported simulations
./set_simulation_parameters.py

# Link to build the merger tree
tangos link

# Import AHF Halo catalogue properties
tangos import-properties `cat import_properties | grep -v '^#'`

# Write galaxy properties
tangos write `cat write_properties | grep -v '^#'`
