#/bin/bash

export TANGOS_SIMULATION_FOLDER=`realpath $1`
export TANGOS_DB_CONNECTION=$2
export TANGOS_PROPERTY_MODULES=properties

# Get my own directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export PYTHONPATH=$SCRIPT_DIR

echo $TANGOS_SIMULATION_FOLDER

# Import the data
for i in $TANGOS_SIMULATION_FOLDER/*
do
    tangos add `basename $i`
done

# Set any non-default simulation properties for all imported simulations
$SCRIPT_DIR/set_simulation_parameters.py

# Link to build the merger tree
tangos link

# Import AHF Halo catalogue properties
tangos import-properties `cat $SCRIPT_DIR/import_properties | grep -v '^#'`

# Write galaxy properties
tangos write `cat $SCRIPT_DIR/write_properties | grep -v '^#'`
