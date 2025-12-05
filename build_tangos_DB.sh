#/bin/bash

if [ "$#" -lt 1 ]
then
    echo 'Please specifiy a config file'
    exit 1
fi

export TANGOS_PROPERTY_MODULES=properties
export $(cat $1 | grep -v '^#' | xargs)
export TANGOS_SIMULATION_FOLDER=$(realpath $TANGOS_SIMULATION_FOLDER)

# Get my own directory
SCRIPT_DIR=$(dirname "$0")

export PYTHONPATH=$SCRIPT_DIR

echo $TANGOS_DB_CONNECTION

if [ -d $TANGOS_SIMULATION_FOLDER ]
then
    echo "Adding Data From" $TANGOS_SIMULATION_FOLDER
else
    echo "Simulation folder does not exist"
    exit 1
fi

# Tangos won't import properties if PYTEST_CURRENT_TEST is set
unset PYTEST_CURRENT_TEST

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
tangos write `cat $SCRIPT_DIR/write_properties | grep -v '^#'` --include-only="NGas()>$MIN_GAS" --include-only="NStar()>$MIN_STAR"
