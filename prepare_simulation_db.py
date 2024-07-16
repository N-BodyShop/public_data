import os
import re
from collections import defaultdict

def is_snapshot_file(filepath, sim_name):
    pattern = f"{sim_name}/{sim_name}.snapshot/{sim_name}.snapshot"
    return filepath.endswith(pattern)

def is_file_readable(filepath, sim_name):
    if not os.path.isfile(filepath):
        return True  # Directories are considered "readable"
    if is_snapshot_file(filepath, sim_name):
        return True  # Snapshot files are exempted from readability check
    try:
        with open(filepath, 'rb') as f:
            f.read(1)
        return True
    except IOError:
        #print(f"Warning: File {filepath} is not readable.")
        return False

def create_symlink(source, dest, sim_name, overwrite=True):
    if not is_file_readable(source, sim_name):
        #print(f"Skipping unreadable file: {source}")
        return

    if os.path.exists(dest):
        if overwrite:
            os.remove(dest)
            os.symlink(source, dest)
            #print(f"Overwritten symlink: {dest} -> {source}")
        #else:
            #print(f"File already exists: {dest}")
    else:
        os.symlink(source, dest)
        #print(f"Created symlink: {dest} -> {source}")


def clean_ahf_filename(filename):
    # Check if it's an AHF file
    if 'AHF' not in filename:
        #print('AHF not in file name')
        return filename

    # Extract the snapshot number
    match = re.search(r'\b(00\d{3,4})\b', filename)
    # print(match.group(0))
    if not match:
        # print('could not find snapshot number')
        return filename

    snapshot_number = match.group(1)
    prefix = filename[:match.start()]

    # Check if it's already in the right format
    if re.search(rf'{snapshot_number}\.z\d\.\d\d\d\.AHF_', filename):
        # print('already in right format')
        return filename

    # If not in the right format, clean it up
    parts = filename.split(snapshot_number)
    if len(parts) > 1:
        suffix = parts[1]
        z_match = re.search(rf'\.z\d\.\d\d\d\.AHF_([^\s]+)', suffix)
        if z_match:
            return f"{parts[0]}{snapshot_number}{z_match.group()}"
    # print('could not clean')
    # If we couldn't clean it up, rais an error
    raise ValueError(f"Could not clean AHF filename: {filename}")
    #return filename

def get_file_category(filename):
    # Define patterns for different file categories
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
            if file.endswith('.amiga.grp'):
                continue
            if file.endswith('.amiga.gtp'):
                continue
            if file.endswith('.amiga.stat'):
                continue

            source_file = os.path.join(root, file)
            cleaned_file = clean_ahf_filename(file)
            category = get_file_category(cleaned_file)

            # Store file info for later processing
            file_categories[category].append((source_file, cleaned_file, dest_subdir))

    # Process files by category
    for category in ['snapshot', 'ahf', 'parameter', 'other']:
        for source_file, cleaned_file, dest_subdir in file_categories[category]:
            dest_file = os.path.join(dest_subdir, cleaned_file)
            create_symlink(source_file, dest_file, sim_name)

def process_simulations(match_string, base_source_dir, base_dest_dir):
    for sim_folder in os.listdir(base_source_dir):
        if match_string in sim_folder:
            source_dir = os.path.join(base_source_dir, sim_folder)
            print(f"Processing simulation: {sim_folder}")
            process_simulation(source_dir, base_dest_dir)

if __name__ == '__main__':
    # location of the destination directory
    base_dest_dir = '/data/REPOSITORY/public_data/'


    #source directory for Marvel volumes
    base_source_dir = '/data/REPOSITORY/dwarf_volumes'


    simulations = ['cptmarvel.cosmo25cmb.4096g5HbwK1BH',
                   'elektra.cosmo25cmb.4096g5HbwK1BH',
                   'rogue.cosmo25cmb.4096g5HbwK1BH',
                   'storm.cosmo25cmb.4096g1HsbBH',
                   'storm.cosmo25cmb.4096g5HbwK1BH']

    for sim in simulations:
        #print(f"Processing Marvel simulation: {sim}")
        process_simulations(sim, base_source_dir, base_dest_dir)

    #source directory for gs14
    base_source_dir = '/data/REPOSITORY/public'

    simulations = ['h239.cosmo50cmb.3072g14HMbwK',
                   'h258.cosmo50cmb.3072g14HMbwK',
                   'h277.cosmo50cmb.3072g14HMbwK']
    for sim in simulations:
        process_simulations(sim, base_source_dir, base_dest_dir)

    print("All simulations processed.")