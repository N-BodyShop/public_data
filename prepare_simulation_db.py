# File: prepare_simulation_db.py

#create a new folders in the public_data directory that prevent tangos from trying to add both
# ahf 200 and ahf 178 files to the database.

import os
import re


def create_symlink(source, dest, overwrite=False):
    if os.path.exists(dest):
        if overwrite:
            os.remove(dest)
            os.symlink(source, dest)
            print(f"Overwritten symlink: {dest} -> {source}")
        else:
            print(f"File already exists: {dest}")
    else:
        os.symlink(source, dest)
        print(f"Created symlink: {dest} -> {source}")


def process_simulation(sim_dir, dest_base):
    sim_name = os.path.basename(sim_dir)
    dest_dir = os.path.join(dest_base, sim_name)

    # Create the destination directory if it doesn't exist
    os.makedirs(dest_dir, exist_ok=True)

    for root, dirs, files in os.walk(sim_dir):
        if  'ahf_200' in dirs:
            dirs.remove('ahf_200')

        # Create relative path
        rel_path = os.path.relpath(root, sim_dir)
        dest_subdir = os.path.join(dest_dir, rel_path)

        # Create subdirectory in destination
        os.makedirs(dest_subdir, exist_ok=True)

        # Create symlinks for files
        for file in files:
            source_file = os.path.join(root, file)
            dest_file = os.path.join(dest_subdir, file)
            create_symlink(source_file, dest_file)


def process_simulations(match_string, base_source_dir, base_dest_dir, exclude_ahf_200=False):
    for sim_folder in os.listdir(base_source_dir):
        if match_string in sim_folder:
            source_dir = os.path.join(base_source_dir, sim_folder)
            print(f"Processing simulation: {sim_folder}")
            process_simulation(source_dir, base_dest_dir, exclude_ahf_200)


if __name__ == '__main__':
    base_source_dir = '/data/REPOSITORY/dwarf_volumes'
    base_dest_dir = '/data/REPOSITORY/public_data/Marvel'

    simulations = ['cptmarvel', 'storm', 'rogue', 'elektra']

    for sim in simulations:
        print(f"Processing {sim} simulations...")
        process_simulations(sim, base_source_dir, base_dest_dir)

    print("All simulations processed.")