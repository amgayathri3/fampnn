import os
import csv
from rfantibody.util.pose import Pose

def group_contiguous_positions(positions):
    if not positions:
        return []
    positions = sorted(positions)
    ranges = []
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            ranges.append((start, prev))
            start = prev = pos
    ranges.append((start, prev))
    return ranges

def invert_ranges(total_len, cdr_ranges):
    cdr_set = set()
    for start, end in cdr_ranges:
        cdr_set.update(range(start, end + 1))
    non_cdrs = []
    start = None
    for i in range(1, total_len + 1):
        if i not in cdr_set:
            if start is None:
                start = i
        elif start is not None:
            non_cdrs.append((start, i - 1))
            start = None
    if start is not None:
        non_cdrs.append((start, total_len))
    return non_cdrs

def format_ranges(ranges, chain, valid_residues=None):
    formatted = []
    for start, end in ranges:
        if valid_residues is not None:
            in_pdb = any((chain, i) in valid_residues for i in range(start, end + 1))
            if not in_pdb:
                continue
        formatted.append(f"{chain}{start}-{end}")
    return formatted

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

def generate_csv_for_fampnn(pdb_dir, out_csv, pdb_key_list=None, invert_mask=True):
    rows = []

    keys = None
    if pdb_key_list:
        with open(pdb_key_list) as f:
            keys = set(line.strip() for line in f if line.strip())

    for fname in sorted(os.listdir(pdb_dir)):
        if not fname.endswith(".pdb"):
            continue
        pdb_path = os.path.join(pdb_dir, fname)
        pdb_id = os.path.splitext(fname)[0]

        if keys and pdb_id not in keys:
            continue

        pose = Pose.from_pdb(pdb_path)
        valid_residues = extract_residues_from_pdb(pdb_path)

        print(f"[DEBUG] {pdb_id} extracted CDRs:", pose.cdr_dict)
        print(f"[DEBUG] {pdb_id} valid residues:", sorted(valid_residues))

        fixed_seq_ranges = []

        for chain, cdr_positions in pose.cdr_dict.items():
            grouped_cdrs = group_contiguous_positions(cdr_positions)
            if not grouped_cdrs:
                continue
            seq_len = max(cdr_positions)
            if invert_mask:
                non_cdr_ranges = invert_ranges(seq_len, grouped_cdrs)
                # skip valid_residues filtering to test for mismatch
                fixed_seq_ranges.extend(format_ranges(non_cdr_ranges, chain, valid_residues=None))
            else:
                fixed_seq_ranges.extend(format_ranges(grouped_cdrs, chain, valid_residues=None))

        fixed_seq_str = ",".join(fixed_seq_ranges)
        rows.append([pdb_id, fixed_seq_str, ""]) 

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "fixed_position", "fixed_sc_position"])
        writer.writerows(rows)

# Entry point
if __name__ == "__main__":
    import argparse
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument("--pdb_dir", type=str, required=True, help="Directory with input PDBs")
        parser.add_argument("--out_csv", type=str, required=True, help="Output FAMPnn fixed_pos_csv file")
        parser.add_argument("--pdb_key_list", type=str, default=None, help="Optional key list of PDBs to include")
        return parser.parse_args()

    args = parse_args()
    generate_csv_for_fampnn(args.pdb_dir, args.out_csv, pdb_key_list=args.pdb_key_list, invert_mask=True)
