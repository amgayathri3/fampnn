#!/usr/bin/env python3
import os
import argparse
import subprocess
import torch

from generate_fampnn_cdr_mask_csv import generate_csv_for_fampnn
from rfantibody.rf2.modules.pose_util import (
    pose_from_remarked,
    pose_to_remarked_pdblines,
    pdblines_to_pdb,
    Pose
)

def get_pdb_keys(pdb_dir):
    return [
        os.path.splitext(f)[0]
        for f in sorted(os.listdir(pdb_dir))
        if f.endswith(".pdb")
    ]

def write_key_list(pdb_keys, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for key in pdb_keys:
            f.write(key + "\n")

def run_fampnn(pdb_dir, key_list_path, fixed_pos_csv, out_dir):
    with open(key_list_path) as f:
        keys = [l.strip() for l in f if l.strip()]
    tmp_list = key_list_path.replace(".txt", "_fullpaths.txt")
    with open(tmp_list, "w") as f:
        for k in keys:
            f.write(f"{k}.pdb\n")

    cmd = [
        "python3", "fampnn/inference/seq_design.py",
        "checkpoint_path=weights/fampnn_0_3.pt",
        f"pdb_dir={pdb_dir}",
        f"pdb_key_list={tmp_list}",
        f"fixed_pos_csv={fixed_pos_csv}",
        f"out_dir={out_dir}"
    ]
    print("Running FAMPnn:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def extract_backbone(pdb_path):
    """
    Read N, CA, C coords and sequence from a PDB.
    Returns: xyz Tensor [L,3,3], seq Tensor[L]
    """
    aa3 = [
        'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
        'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
    ]
    aa2idx = {r:i for i,r in enumerate(aa3)}

    coords = []
    seq = []
    last = (None, None)
    current = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resnum = int(line[22:26])
            key = (chain, resnum)
            if key != last:
                if current:
                    coords.append(current)
                    seq.append(last_resn)
                current = []
                last = key
            atom_name = line[12:16].strip()
            x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
            if atom_name == "N":
                current.insert(0, (x, y, z))
            elif atom_name == "CA":
                current.insert(1, (x, y, z))
                last_resn = line[17:20].strip()
            elif atom_name == "C":
                current.insert(2, (x, y, z))
        if current:
            coords.append(current)
            seq.append(last_resn)

    L = len(coords)
    xyz = torch.zeros((L, 3, 3), dtype=torch.float32)
    s = torch.zeros((L,), dtype=torch.long)
    for i, atom_list in enumerate(coords):
        for j in range(3):
            xyz[i, j] = torch.tensor(atom_list[j])
        s[i] = aa2idx.get(seq[i], 0)
    return xyz, s

def restore_chain_ids(pdb_dir, out_dir):
    for fname in sorted(os.listdir(out_dir)):
        if not fname.endswith(".pdb") or "_sample" not in fname:
            continue

        sample_path = os.path.join(out_dir, fname)
        tag = fname.split("_sample")[0]
        orig_path = os.path.join(pdb_dir, f"{tag}.pdb")
        if not os.path.exists(orig_path):
            print(f"Warning: original PDB not found for {fname}")
            continue

        try:
            # load original pose (with H/L/T and CDR masks)
            pose = pose_from_remarked(orig_path)

            # extract new backbone coords and sequence
            new_xyz, new_seq = extract_backbone(sample_path)

            # rebuild pose with new data
            restored = Pose(
                xyz=new_xyz,
                seq=new_seq,
                atom_mask=pose.atom_mask,
                cdrs=pose.cdrs,
                idx=pose.idx,
                chain_dict=pose.chain_dict
            )

            # regenerate HLT PDB with REMARKs
            lines = pose_to_remarked_pdblines(restored)
            pdblines_to_pdb(lines, sample_path)

            print(f"Restored H/L/T chains for {fname}")
        except Exception as e:
            print(f"Error restoring {fname}: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", required=True, help="Directory of input HLT PDBs")
    parser.add_argument("--out_dir", required=True, help="Directory to save FAMPnn outputs")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    key_list = os.path.join("examples/pdb_key_lists", f"{os.path.basename(args.pdb_dir)}.txt")
    csv_file = os.path.join("examples/fixed_pos_csvs", f"{os.path.basename(args.pdb_dir)}_not_cdr_mask.csv")

    write_key_list(get_pdb_keys(args.pdb_dir), key_list)
    generate_csv_for_fampnn(args.pdb_dir, csv_file)

    print("Starting FAMPnn design...")
    run_fampnn(args.pdb_dir, key_list, csv_file, args.out_dir)

    print("Post-processing outputs to HLT format...")
    restore_chain_ids(args.pdb_dir, args.out_dir)

if __name__ == "__main__":
    main()

