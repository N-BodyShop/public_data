#!/bin/bash

export TANGOS_SIMULATION_FOLDER="/home/bk639/data/MerianSIDM"
export TANGOS_DB_CONNECTION="${TANGOS_SIMULATION_FOLDER}/MerianSIDM.db"

echo "Using database: $TANGOS_DB_CONNECTION"

# Loop through immediate subfolders in the directory
for folder in "${TANGOS_SIMULATION_FOLDER}"/*/; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        echo "Processing folder: $folder_name"

        if [ "$folder_name" = "adiabatic" ]; then
            # For the adiabatic folder, process its subfolders
            for sim in "$folder"*/; do
                if [ -d "$sim" ]; then
                    sim_name=$(basename "$sim")
                    echo "Processing adiabatic simulation: $sim_name"

                    # Add the simulation using adiabatic/simname format
                    tangos add "adiabatic/${sim_name}" --backend multiprocessing-32 --no-renumber

                    if [ $? -eq 0 ]; then
                        echo "Successfully added adiabatic/${sim_name} to the database"
                    else
                        echo "Error: Failed to add adiabatic/${sim_name} to the database"
                    fi
                fi
            done
        else
            # Process other folders normally
            tangos add "${folder_name}" --backend multiprocessing-32 --no-renumber

            if [ $? -eq 0 ]; then
                echo "Successfully added $folder_name to the database"
            else
                echo "Error: Failed to add $folder_name to the database"
            fi
        fi
    fi
done

echo "All simulations processed"