#!/bin/bash

#this script must be run before write_import.sh

# Base directory for Marvel simulations
export TANGOS_SIMULATION_FOLDER="/data/REPOSITORY/public_data"

# Python path export (uncomment if needed)
# export PYTHONPATH="/home/bk639/:$PYTHONPATH"

# Function to create database connection string
create_db_connection() {
    local sim_name=$1
    echo "${TANGOS_SIMULATION_FOLDER}/${sim_name}.db"
}

# Properties to import from import_properties.txt file
# Read properties from file into a variable, removing any new lines
properties=$(tr '\n' ' ' < import_properties | sed 's/ $//')

echo "Properties to import: $properties"


# Loop through immediate subfolders in the Marvel directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Processing folder: $folder_name"

        # load existing database
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION"

       # Construct and execute the tangos command with properties
        tangos_command="tangos import-properties $properties --sim \"$folder_name\" --backend multiprocessing-4 "
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
