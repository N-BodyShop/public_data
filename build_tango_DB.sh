#/usr/bin/env bash

# Import the data
#for i in $TANGOS_SIMULATION_FOLDER/MUGS2*
#do
#    tangos add `basename $i`
#done
tangos add MUGS2_g1536

# Link to build the merger tree
tangos link

# Import AHF Halo catalogue properties
tangos import-properties `cat import_properties`
