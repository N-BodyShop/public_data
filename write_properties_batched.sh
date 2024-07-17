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

# Define property batches
declare -a batches=(
    "shrink_center finder_mass finder_dm_mass finder_gas_mass finder_star_mass" #instantaneous_SFR"
    "dm_density_profile dm_mass_profile gas_density_profile gas_mass_profile star_density_profile star_mass_profile"# tot_mass_profile"
    "star_metal_profile gas_metal_profile star_Fe_profile gas_Fe_profile star_Ox_profile gas_Ox_profile cold_gas_metal_profile cold_gas_Fe_profile cold_gas_Ox_profile"
    "cold_gas_mass_profile warm_gas_mass_profile hot_gas_mass_profile"
    "vrdisp_stars vrdisp_gas vrdisp_dm vrdisp_stars_3d vrdisp_gas_3d vrdisp_dm_3d vrdisp_encl_stars vrdisp_encl_gas vrdisp_encl_dm vrdisp_encl_stars_3d vrdisp_encl_gas_3d vrdisp_encl_dm_3d"
    "v_surface_brightness b_surface_brightness i_surface_brightness"
    "SFR_histogram"# outflow_rates"
)

# Loop through immediate subfolders in the sim directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Processing folder: $folder_name"

        # Load existing database
        export TANGOS_DB_CONNECTION=$(create_db_connection "$folder_name")
        echo "Using database: $TANGOS_DB_CONNECTION"

        # Process each batch
        for batch in "${batches[@]}"; do
            echo "Processing batch: $batch"
            
            # Construct and execute the tangos write command
            tangos_command="tangos write $batch --sim \"$folder_name\" --backend multiprocessing-32 --load-mode=server-shared-mem"
            echo "Executing: $tangos_command"
            eval $tangos_command

            # Check if the command was successful
            if [ $? -eq 0 ]; then
                echo "Successfully added properties for $folder_name to the database"
            else
                echo "Error: Failed to add properties for $folder_name to the database"
                break  # Exit the batch loop if there's an error
            fi
        done
    fi
done

echo "All simulations processed"
