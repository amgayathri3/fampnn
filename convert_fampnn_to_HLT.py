import argparse
import csv
from rfantibody.util.pose import Pose

def parse_args():
    parser = argparse.ArgumentParser(description="Convert PDB to HLT format and annotate CDRs using FAMPnn mask CSV")
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("--heavy", "-H", required=True, help="Heavy chain ID in input")
    parser.add_argument("--light", "-L", required=True, help="Light chain ID in input")
    parser.add_argument("--target", "-T", required=True, help="Target chain ID(s), comma-separated")
    parser.add_argument("--output", "-o", help="Output PDB file (HLT format)")
    parser.add_argument("--mask_csv", required=True, help="Path to FAMPnn fixed_pos_csv file")
    return parser.parse_args()

def remap_chains(pose, chain_map):
    for i, old_chain in enumerate(pose.chain):
        if old_chain in chain_map:
            pose.chain[i] = chain_map[old_chain]
    if hasattr(pose, 'res_chain'):
        for i, old_chain in enumerate(pose.res_chain):
            if old_chain in chain_map:
                pose.res_chain[i] = chain_map[old_chain]

def reorder_chains(pose, desired_order):
    order_dict = {chain: i for i, chain in enumerate(desired_order)}
    indices = sorted(range(len(pose.chain)), key=lambda i: order_dict.get(pose.chain[i], 999))
    for attr in ['chain', 'resname', 'resseq', 'icode', 'coords', 'atomname', 'b', 'q']:
        if hasattr(pose, attr):
            arr = getattr(pose, attr)
            if isinstance(arr, list):
                setattr(pose, attr, [arr[i] for i in indices])
            else:  # numpy array
                setattr(pose, attr, arr[indices])

def parse_fixed_positions(fixed_str):
    positions = set()
    for part in fixed_str.split(','):
        if not part:
            continue
        chain = part[0]
        rng = part[1:]
        if '-' in rng:
            start, end = map(int, rng.split('-'))
            for i in range(start, end + 1):
                positions.add((chain, i))
        else:
            positions.add((chain, int(rng)))
    return positions

def extract_residues_from_pdb(pdb_path):
    residues = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and len(line) >= 26:
                chain = line[21]
                try:
                    resnum = int(line[22:26])
                    residues.add((chain, resnum))
                except ValueError:
                    continue
    return residues

def get_cdr_residues_from_csv(csv_path, pdb_name, all_residues):
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['name'] == pdb_name:
                fixed = parse_fixed_positions(row['fixed_position'])
                return all_residues - fixed
    return set()

def main():
    args = parse_args()
    target_chains = [x.strip() for x in args.target.split(',')]
    out_pdb = args.output or args.input_pdb.replace('.pdb', '_HLT.pdb')

    # 1. Load PDB as Pose
    pose = Pose.from_pdb(args.input_pdb)

    # 2. Remap chain IDs
    chain_map = {args.heavy: 'H', args.light: 'L'}
    for t in target_chains:
        chain_map[t] = 'T'
    remap_chains(pose, chain_map)

    # 3. Reorder chains: H, L, T
    reorder_chains(pose, ['H', 'L', 'T'])

    # 4. Write out PDB
    pose.dump_pdb(out_pdb)

    # 5. Annotate CDRs with REMARKs using FAMPnn mask CSV
    pdb_name = out_pdb.split('/')[-1].replace('_HLT.pdb', '').replace('.pdb', '')
    all_residues = extract_residues_from_pdb(out_pdb)
    cdr_residues = get_cdr_residues_from_csv(args.mask_csv, pdb_name, all_residues)

    with open(out_pdb, 'a') as f:
        for chain, resnum in sorted(cdr_residues):
            f.write(f"REMARK CDR CHAIN {chain} RES {resnum}\n")

if __name__ == '__main__':
    main()
