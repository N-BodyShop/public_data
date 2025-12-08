#/bin/bash

if [ "$#" -lt 1 ]
then
    echo 'Please specifiy a config file'
    exit 1
fi

# Load and process the config file
export TANGOS_PROPERTY_MODULES=properties
export $(cat $1 | grep -v '^#' | xargs)
export TANGOS_SIMULATION_FOLDER=$(realpath $TANGOS_SIMULATION_FOLDER)


if [ -z $MIN_GAS ]
then
    MIN_GAS=0
fi

if [ -z $MIN_STAR ]
then
    MIN_STAR=0
fi

if [ -z $NPROCS ]
then
    echo "Running in serial mode.  This may be slow!"
else
    if [ "$MPI" = true ]
    then
        PARALLEL="--backend=mpi4py"
        RUNNER="mpirun -np "$NPROCS
        if [ "$SERVER" = true ]
        then
            LOAD_MODE="--load-mode=server"
        fi
    else
        PARALLEL="--backend=multiprocessing-"$NPROCS
        if [ "$SERVER" = true ]
        then
            LOAD_MODE="--load-mode=server-shared-mem"
        fi
    fi
fi


# Get my own directory
SCRIPT_DIR=$(dirname "$0")

export PYTHONPATH=$SCRIPT_DIR

echo $TANGOS_DB_CONNECTION

if [ -d $TANGOS_SIMULATION_FOLDER ]
then
    echo "Adding Data From" $TANGOS_SIMULATION_FOLDER
else
    echo "Simulation folder does not exist"
    exit 2
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
$RUNNER tangos link $PARALLEL

# Import AHF Halo catalogue properties
tangos import-properties `cat $SCRIPT_DIR/import_properties | grep -v '^#'`

# Write galaxy properties
$RUNNER tangos write `cat $SCRIPT_DIR/write_properties | grep -v '^#'` --include-only="NGas()>$MIN_GAS" --include-only="NStar()>$MIN_STAR" $PARALLEL $LOAD_MODE

