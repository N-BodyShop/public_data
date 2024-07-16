#!/bin/bash

# Base directory for Marvel simulations
export TANGOS_SIMULATION_FOLDER="/data/REPOSITORY/public_data/Marvel"

# Python path export (uncomment if needed)
# export PYTHONPATH="/home/bk639/:$PYTHONPATH"

# Function to create database connection string
create_db_connection() {
    local sim_name=$1
    echo "${TANGOS_SIMULATION_FOLDER}/${sim_name}.db"
}


# Loop through immediate subfolders in the Marvel directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Processing folder: $folder_name"

        # Create a unique database file for this simulation in the Marvel directory
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION"

        # Add the simulation to the database using the folder name relative to TANGOS_SIMULATION_FOLDER
        tangos add "${folder_name}"

        # Check if the command was successful
        if [ $? -eq 0 ]; then
            echo "Successfully added $folder_name to the database"

            #import basic properties
            echo "Importing properties for $folder_name"
            tangos import-properties --sim "${folder_name}"
            echo "Properties imported for $folder_name"
            # Create links for this simulation
            echo "Creating links for $folder_name"
            tangos link --sim "${folder_name}"
            echo "Links created for $folder_name"
        else
            echo "Error: Failed to add $folder_name to the database"
        fi

    fi
done

echo "All simulations processed"