#!/bin/bash

# Base directory for Marvel simulations
export TANGOS_SIMULATION_FOLDER="/data/REPOSITORY/public_data"

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

        # Set the database connection
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION"

        # Add the simulation to the database using the folder name relative to TANGOS_SIMULATION_FOLDER
        tangos link "${folder_name}" --backend multiprocessing-32

        # Check if the command was successful
        if [ $? -eq 0 ]; then
            echo "Successfully added $folder_name to the database"

        else
            echo "Error: Failed to add $folder_name to the database"
        fi

    fi
done

echo "All simulations processed"