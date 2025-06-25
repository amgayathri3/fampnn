#!/usr/bin/env python3
import os
import argparse
import subprocess
from generate_fampnn_cdr_mask_csv import generate_csv_for_fampnn

def get_pdb_keys(pdb_dir):
    keys = []
    for fname in sorted(os.listdir(pdb_dir)):
        if fname.endswith(".pdb"):
            keys.append(os.path.splitext(fname)[0])
    return keys

def write_key_list(pdb_keys, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for key in pdb_keys:
            f.write(f"{key}\n")

def run_fampnn(pdb_dir, key_list_path, fixed_pos_csv_path, out_dir, temperature=1.0, num_samples=1):
    with open(key_list_path) as f:
        keys = [line.strip() for line in f if line.strip()]
    tmp_key_list = key_list_path.replace(".txt", "_fullpaths.txt")
    with open(tmp_key_list, "w") as f:
        for k in keys:
            f.write(f"{k}.pdb\n")

    cmd = [
        "python3", os.path.abspath("fampnn/inference/seq_design.py"),
        f"checkpoint_path={os.path.abspath('weights/fampnn_0_3.pt')}",
        f"pdb_dir={os.path.abspath(pdb_dir)}",
        f"pdb_key_list={os.path.abspath(tmp_key_list)}",
        f"fixed_pos_csv={os.path.abspath(fixed_pos_csv_path)}",
        f"out_dir={os.path.abspath(out_dir)}",
        f"temperature={temperature}",
        f"num_seqs_per_pdb={num_samples}"
    ]

    print("[INFO] Running FAMPnn command:")
    print(" ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] FAMPnn failed with: {e}")
        exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory of input PDB files")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save FAMPnn outputs")
    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature (default: 1.0)")
    parser.add_argument("--num_samples", type=int, default=1, help="Number of samples per input structure (default: 1)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pdb_keys = get_pdb_keys(args.pdb_dir)
    key_list_path = os.path.join("examples/pdb_key_lists", f"{os.path.basename(args.pdb_dir)}_key_list.txt")
    write_key_list(pdb_keys, key_list_path)
    print(f"[INFO] Saved key list to {key_list_path}")

    fixed_pos_csv_path = os.path.join("examples/fixed_pos_csvs", f"{os.path.basename(args.pdb_dir)}_not_cdr_mask.csv")
    os.makedirs(os.path.dirname(fixed_pos_csv_path), exist_ok=True)
    generate_csv_for_fampnn(
        pdb_dir=args.pdb_dir,
        out_csv=fixed_pos_csv_path
    )
    print(f"[INFO] Saved fixed_pos_csv to {fixed_pos_csv_path}")

    print("[INFO] Running FAMPnn sequence design...")
    run_fampnn(
        pdb_dir=args.pdb_dir,
        key_list_path=key_list_path,
        fixed_pos_csv_path=fixed_pos_csv_path,
        out_dir=args.out_dir,
        temperature=args.temperature,
        num_samples=args.num_samples
    )

if __name__ == "__main__":
    main()

