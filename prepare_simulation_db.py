import os
import re
from collections import defaultdict

'''
This script creates symlinks for files in simulation directories.
Key features:
1. Skips ahf_200 directories and files ending in .amiga.grp, .amiga.gtp, and .amiga.stat.
2. Avoids unreadable files.
3. Cleans AHF filenames to the format snapshot_number.zX.YYY.AHF_* for tangos database compatibility.
4. Handles various directory structures (root or snapshot subdirectory).
5. Overwrites duplicate files by default (configurable).
6. Processes both Marvel and gs14 simulations.
'''

def is_snapshot_file(filepath, sim_name):
    # Identifies snapshot files based on their path pattern
    pattern = f"{sim_name}/{sim_name}.snapshot/{sim_name}.snapshot"
    return filepath.endswith(pattern)

def is_file_readable(filepath, sim_name):
    # Checks if a file is readable, considering directories and snapshot files as always readable
    if not os.path.isfile(filepath):
        return True  # Directories are considered "readable"
    if is_snapshot_file(filepath, sim_name):
        return True  # Snapshot files are exempted from readability check
    try:
        with open(filepath, 'rb') as f:
            f.read(1)
        return True
    except IOError:
        return False

def create_symlink(source, dest, sim_name, overwrite=True):
    # Creates a symlink, optionally overwriting existing files
    if not is_file_readable(source, sim_name):
        return

    if os.path.exists(dest):
        if overwrite:
            os.remove(dest)
            os.symlink(source, dest)
    else:
        os.symlink(source, dest)

def clean_ahf_filename(filename):
    # Cleans AHF filenames to the format required by the tangos database
    if 'AHF' not in filename:
        return filename

    match = re.search(r'\b(00\d{3,4})\b', filename)
    if not match:
        return filename

    snapshot_number = match.group(1)
    prefix = filename[:match.start()]

    if re.search(rf'{snapshot_number}\.z\d\.\d\d\d\.AHF_', filename):
        return filename

    parts = filename.split(snapshot_number)
    if len(parts) > 1:
        suffix = parts[1]
        z_match = re.search(rf'\.z\d\.\d\d\d\.AHF_([^\s]+)', suffix)
        if z_match:
            return f"{parts[0]}{snapshot_number}{z_match.group()}"

    raise ValueError(f"Could not clean AHF filename: {filename}")

def get_file_category(filename):
    # Categorizes files based on their names
    patterns = {
        'snapshot': r'\.snapshot$',
        'ahf': r'AHF',
        'amiga': r'\.amiga$',
        'parameter': r'\.(param|txt)$',
        'other': r'.*'
    }

    for category, pattern in patterns.items():
        if re.search(pattern, filename):
            return category
    return 'other'

def process_simulation(sim_dir, dest_base):
    # Processes a single simulation directory, creating symlinks for all relevant files
    sim_name = os.path.basename(sim_dir)
    dest_dir = os.path.join(dest_base, sim_name)
    os.makedirs(dest_dir, exist_ok=True)

    file_categories = defaultdict(list)

    for root, dirs, files in os.walk(sim_dir):
        if 'ahf_200' in dirs:
            dirs.remove('ahf_200')

        rel_path = os.path.relpath(root, sim_dir)
        dest_subdir = os.path.join(dest_dir, rel_path)
        os.makedirs(dest_subdir, exist_ok=True)

        for file in sorted(files):
            if file.endswith(('.amiga.grp', '.amiga.gtp', '.amiga.stat')):
                continue

            source_file = os.path.join(root, file)
            cleaned_file = clean_ahf_filename(file)
            category = get_file_category(cleaned_file)

            file_categories[category].append((source_file, cleaned_file, dest_subdir))

    # Process files by category
    for category in ['snapshot', 'ahf', 'parameter', 'other']:
        for source_file, cleaned_file, dest_subdir in file_categories[category]:
            dest_file = os.path.join(dest_subdir, cleaned_file)
            create_symlink(source_file, dest_file, sim_name)

def process_simulations(match_string, base_source_dir, base_dest_dir):
    # Processes all simulation directories that match the given string
    for sim_folder in os.listdir(base_source_dir):
        if match_string in sim_folder:
            source_dir = os.path.join(base_source_dir, sim_folder)
            print(f"Processing simulation: {sim_folder}")
            process_simulation(source_dir, base_dest_dir)

if __name__ == '__main__':
    # Destination directory for all processed simulations
    base_dest_dir = '/data/REPOSITORY/public_data/'

    # Process Marvel volumes
    base_source_dir = '/data/REPOSITORY/dwarf_volumes'
    marvel_simulations = [
        'cptmarvel.cosmo25cmb.4096g5HbwK1BH',
        'elektra.cosmo25cmb.4096g5HbwK1BH',
        'rogue.cosmo25cmb.4096g5HbwK1BH',
        'storm.cosmo25cmb.4096g1HsbBH',
        'storm.cosmo25cmb.4096g5HbwK1BH'
    ]

    for sim in marvel_simulations:
        process_simulations(sim, base_source_dir, base_dest_dir)

    # Process gs14 simulations
    base_source_dir = '/data/REPOSITORY/public'
    gs14_simulations = [
        'h239.cosmo50cmb.3072g14HMbwK',
        'h258.cosmo50cmb.3072g14HMbwK',
        'h277.cosmo50cmb.3072g14HMbwK'
    ]

    for sim in gs14_simulations:
        process_simulations(sim, base_source_dir, base_dest_dir)

    print("All simulations processed.")