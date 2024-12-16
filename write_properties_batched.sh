#!/bin/bash

# Set environment variables
export TANGOS_SIMULATION_FOLDER="/home/bk639/data/public_test"
export PYTHONPATH="/home/bk639/public_data/"
export TANGOS_PROPERTY_MODULES=properties
export TANGOS_DB_CONNECTION="/home/bk639/data/public_test/test.db"

# Function to create database connection string
create_db_connection() {
    local sim_name=$1
    echo "${TANGOS_SIMULATION_FOLDER}/${sim_name}.db"
}

# Function to log messages
log_message() {
    local sim_name=$1
    local batch_name=$2
    local message=$3
    local log_folder="${TANGOS_SIMULATION_FOLDER}/${sim_name}_logs"
    local log_file="${log_folder}/${batch_name}_log.txt"

    # Create log folder if it doesn't exist
    mkdir -p "$log_folder"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$log_file"
    echo "$message"
}

# Define property batches
declare -a write_functions=(
    "shrink_center"
    "finder_mass"
    "finder_dm_mass"
    "SFR_histogram"
    "dm_density_profile"
    "gas_density_profile"
    "star_metal_profile"
    "cold_gas_metal_profile"
    "cold_gas_mass_profile"
    "vrdisp_stars"
    "vrdisp_encl_stars"
    "radius"
    "v_surface_brightness"
    "half_light"
)

# Loop through immediate subfolders in the sim directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        log_message "$folder_name" "main" "Processing folder: $folder_name"

        # Load existing database

        #log_message "$folder_name" "main" "Using database: $TANGOS_DB_CONNECTION"

        # iterate over each function in write_functions
        for function in "${write_functions[@]}"; do
            batch_name=$(echo "$function" | cut -d' ' -f1)  # Use the first function as the batch name
            log_message "$folder_name" "$batch_name" "Processing batch: $function"

            # Construct and execute the tangos write command
            tangos_command="tangos write $function --latest --sim \"$folder_name\" --backend multiprocessing-5 --load-mode=server-shared-mem"
            log_message "$folder_name" "$batch_name" "Executing: $tangos_command"

            # Execute the command and capture its output and exit status
            output=$(eval $tangos_command 2>&1)
            exit_status=$?

            # Log the command output
            log_message "$folder_name" "$batch_name" "Command output: $output"

            # Check if the command was successful
            if [ $exit_status -eq 0 ]; then
                log_message "$folder_name" "$batch_name" "Successfully added properties for $folder_name to the database"
            else
                log_message "$folder_name" "$batch_name" "Error: Failed to add properties for $folder_name to the database"
                break  # Exit the batch loop if there's an error
            fi
        done
    fi
done

log_message "all" "main" "All simulations processed"