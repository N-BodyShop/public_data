import os
import re
from collections import defaultdict


def is_snapshot_file(filepath, sim_name):
    pattern = f"{sim_name}/{sim_name}.snapshot/{sim_name}.snapshot"
    return filepath.endswith(pattern)


def is_file_readable(filepath, sim_name):
    if not os.path.isfile(filepath):
        return True
    if is_snapshot_file(filepath, sim_name):
        return True
    try:
        with open(filepath, 'rb') as f:
            f.read(1)
        return True
    except IOError:
        return False


def create_symlink(source, dest, sim_name, overwrite=True):
    if not is_file_readable(source, sim_name):
        return

    if os.path.exists(dest):
        if overwrite:
            os.remove(dest)
            os.symlink(source, dest)
    else:
        os.symlink(source, dest)


def clean_ahf_filename(filename):
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

    return filename


def get_file_category(filename):
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


def should_skip_file(filename):
    # Skip files ending with .1.iord pattern (duplicates)
    #if re.search(r'\.\d+\.(iord|igasorder|age|massform)$', filename):
     #   return True

    # Skip specific amiga files
    # if filename.endswith(('.amiga.grp', '.amiga.gtp', '.amiga.stat')):
    #     return True

    

    return False


def get_snapshot_number(filename):
    # Extract snapshot number from filename
    match = re.search(r'\b(00\d{3,4})\b', filename)
    return match.group(1) if match else None


def process_simulation(sim_dir, dest_base):
    sim_name = os.path.basename(sim_dir)
    dest_dir = os.path.join(dest_base, sim_name)
    os.makedirs(dest_dir, exist_ok=True)

    file_categories = defaultdict(list)
    processed_snapshots = set()  # Track which snapshot numbers we've processed

    for root, dirs, files in os.walk(sim_dir):
        if 'ahf_200' in dirs:
            dirs.remove('ahf_200')

        for file in sorted(files):
            if should_skip_file(file):
                continue

            source_file = os.path.join(root, file)
            cleaned_file = clean_ahf_filename(file)

            # Get snapshot number if present
            snapshot_num = get_snapshot_number(cleaned_file)

            # If this is a snapshot-related file, check if we've already processed it
            if snapshot_num:
                snapshot_key = f"{snapshot_num}_{os.path.splitext(cleaned_file)[1]}"
                if snapshot_key in processed_snapshots:
                    continue
                processed_snapshots.add(snapshot_key)

            category = get_file_category(cleaned_file)
            file_categories[category].append((source_file, cleaned_file))

    # Process files by category
    for category in ['snapshot', 'ahf', 'parameter', 'other']:
        for source_file, final_name in file_categories[category]:
            dest_file = os.path.join(dest_dir, final_name)
            create_symlink(source_file, dest_file, sim_name)


def process_simulations(match_string, base_source_dir, base_dest_dir):
    for sim_folder in os.listdir(base_source_dir):
        if match_string in sim_folder:
            source_dir = os.path.join(base_source_dir, sim_folder)
            print(f"Processing simulation: {sim_folder}")
            process_simulation(source_dir, base_dest_dir)


if __name__ == '__main__':
    # Marvel simulations
    base_source_dir = '/data/REPOSITORY/dwarf_volumes'
    base_dest_dir = '/home/bk639/data/public_test/Marvel'
    #base_dest_dir = '/home/ae589/Tangos_DB'
    marvel_simulations = [
        'cptmarvel.cosmo25cmb.4096g5HbwK1BH',
        'elektra.cosmo25cmb.4096g5HbwK1BH',
        'rogue.cosmo25cmb.4096g5HbwK1BH',
        #'storm.cosmo25cmb.4096g1HsbBH',
        'storm.cosmo25cmb.4096g5HbwK1BH'


        #'storm.cosmo25cmb.4096'
        

    ]

    for sim in marvel_simulations:
        process_simulations(sim, base_source_dir, base_dest_dir)

    # DCJL simulations
    base_source_dir = '/data/REPOSITORY/e12Gals'
    base_dest_dir = '/home/bk639/data/public_test/DCJL'
    gs14_simulations = [
        'h239.cosmo50cmb.3072g14HMbwK',
        'h258.cosmo50cmb.3072g14HMbwK',
        'h277.cosmo50cmb.3072g14HMbwK'
    ]

    for sim in gs14_simulations:
        process_simulations(sim, base_source_dir, base_dest_dir)

    print("All simulations processed.")