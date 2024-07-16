#!/bin/bash

#this script must be run before write_import.sh

# Base directory for Marvel simulations
export TANGOS_SIMULATION_FOLDER="/data/REPOSITORY/public_data"

export TANGOS_PROPERTY_MODULES=properties

# Python path export (uncomment if needed)
# export PYTHONPATH="/home/bk639/:$PYTHONPATH"

# Function to create database connection string
create_db_connection() {
    local sim_name=$1
    echo "${TANGOS_SIMULATION_FOLDER}/${sim_name}.db"
}

# Read functions from file, remove in-progress ones, and store in a variable
write_functions=$(grep -v '\*\*\*$' write_properties | tr '\n' ' ' | sed 's/ $//')
echo  "Writing properties: $write_functions"


# Loop through immediate subfolders in the Marvel directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Processing folder: $folder_name"

        # load existing database
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION"

       # Construct and execute the tangos write command
        tangos_command="tangos write $write_functions --sim \"$folder_name\" --backend multiprocessing-32 --load-mode=server-shared-mem"
        echo "Executing: $tangos_command"
        eval $tangos_command

        # Check if the command was successful
        if [ $? -eq 0 ]; then
            echo "Successfully added properties for $folder_name to the database"

        else
            echo "Error: Failed to add properties for $folder_name to the database"
        fi

    fi
done

echo "All simulations processed"