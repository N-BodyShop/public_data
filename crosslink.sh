#!/bin/bash

# Array of simulations to process
simulations=(
    "r634.romulus25si2s50v35"
    "r492.romulus25si2s50v35"
    "r468.romulus25si2s50v35"
    "r488.romulus25si2s50v35"
    "r544.romulus25si2s50v35"
    "r597.romulus25si2s50v35"
    "r523.romulus25si2s50v35"
    "r618.romulus25si2s50v35"
    "r431.romulus25cvdXsec"
    "r492.romulus25cvdXsec"
)

# Function to extract the simulation number
get_sim_number() {
    echo "$1" | cut -d'.' -f1
}

# Loop through each simulation and perform the crosslink
for sim in "${simulations[@]}"; do
    sim_number=$(get_sim_number "$sim")
    plain_sim="${sim_number}.romulus25.3072g1HsbBH"
    
    echo "Processing: $sim"
    tangos crosslink "$plain_sim" "$sim" --backend multiprocessing-10
    echo "Completed: $sim"
    echo
done

echo "All crosslinking tasks completed."
