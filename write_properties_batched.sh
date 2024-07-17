#!/bin/bash

# Set environment variables
export TANGOS_SIMULATION_FOLDER="/data/REPOSITORY/public_data"
export PYTHONPATH="/home/bk639/public_data/"
export TANGOS_PROPERTY_MODULES=properties

# Function to create database connection string
create_db_connection() {
    local sim_name=$1
    echo "${TANGOS_SIMULATION_FOLDER}/${sim_name}.db"
}

# Create a logs directory if it doesn't exist
mkdir -p logs

# Read in write functions from write_properties and store as an iterable variable
write_functions=$(grep -v '\*\*\*$' write_properties | tr '\n' ' ' | sed 's/ $//')

# Loop through immediate subfolders in the sim directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Processing folder: $folder_name"

        # Create a directory for this simulation's logs
        sim_log_dir="logs/${folder_name}"
        mkdir -p "$sim_log_dir"

        # Create a main log file for this simulation
        main_log_file="${sim_log_dir}/main_log.txt"
        echo "Processing folder: $folder_name" > "$main_log_file"

        # Load existing database
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION" | tee -a "$main_log_file"

        # iterate over each function in write_functions
        for function in $write_functions; do
            echo "Processing function: $function" | tee -a "$main_log_file"

            # Create a log file for this function
            function_log_file="${sim_log_dir}/${function}_log.txt"
            echo "Processing function: $function for simulation: $folder_name" > "$function_log_file"

            # Construct and execute the tangos write command
            #remove the --latest flag to write all properties to all timesteps
            tangos_command="tangos write $function --latest --sim \"$folder_name\" --backend multiprocessing-32 --load-mode=server-shared-mem"
            echo "Executing: $tangos_command" | tee -a "$main_log_file" "$function_log_file"
            eval $tangos_command 2>&1 | tee -a "$main_log_file" "$function_log_file"

            # Check if the command was successful
            if [ $? -eq 0 ]; then
                echo "Added properties for $folder_name to the database with no critical errors" | tee -a "$main_log_file" "$function_log_file"
            else
                echo "Error: Failed to add properties for $folder_name to the database" | tee -a "$main_log_file" "$function_log_file"
            fi

            echo "----------------------------------------" | tee -a "$main_log_file" "$function_log_file"
        done
    fi
done

echo "All simulations processed"