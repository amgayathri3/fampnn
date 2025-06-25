import os
import argparse
import subprocess
from generate_fampnn_cdr_mask_csv import generate_csv_for_fampnn
from rfantibody.util.pose import Pose
from restore_chain_ids_from_pose import restore_chains

def get_pdb_keys(pdb_dir):
    keys = []
    for fname in sorted(os.listdir(pdb_dir)):
        if fname.endswith(".pdb"):
            keys.append(os.path.splitext(fname)[0])
    return keys

def write_key_list(pdb_keys, output_path):
    with open(output_path, "w") as f:
        for key in pdb_keys:
            f.write(f"{key}\n")

def run_fampnn(pdb_dir, key_list_path, fixed_pos_csv_path, out_dir):
    with open(key_list_path) as f:
        keys = [line.strip() for line in f if line.strip()]
    full_filenames = [f"{k}.pdb" for k in keys]

    tmp_key_list = key_list_path.replace(".txt", "_fullpaths.txt")
    with open(tmp_key_list, "w") as f:
        for fname in full_filenames:
            f.write(fname + "\n")

    cmd = [
        "python3", "fampnn/inference/seq_design.py",
        f"checkpoint_path=weights/fampnn_0_3.pt",
        f"pdb_dir={pdb_dir}",
        f"pdb_key_list={tmp_key_list}",
        f"fixed_pos_csv={fixed_pos_csv_path}",
        f"out_dir={out_dir}"
    ]
    subprocess.run(cmd, check=True)

def restore_chain_ids(pdb_dir, out_dir):
    for fname in os.listdir(out_dir):
        if not fname.endswith(".pdb") or "_sample" not in fname:
            continue
        sample_path = os.path.join(out_dir, fname)
        pdb_id = fname.split("_sample")[0]
        input_pdb = os.path.join(pdb_dir, f"{pdb_id}.pdb")

        if not os.path.exists(input_pdb):
            print(f"[WARN] Skipping {fname}: input PDB not found")
            continue

        try:
            restore_chains(input_pdb, sample_path, sample_path)
            print(f"[INFO] Restored chain IDs: {fname}")
        except Exception as e:
            print(f"[ERROR] Chain restoration failed for {fname}: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory of input PDB files")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save FAMPnn outputs")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pdb_keys = get_pdb_keys(args.pdb_dir)
    key_list_path = os.path.join("examples/pdb_key_lists", f"{os.path.basename(args.pdb_dir)}_key_list.txt")
    os.makedirs(os.path.dirname(key_list_path), exist_ok=True)
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
        out_dir=args.out_dir
    )

    print("[INFO] Restoring original chain IDs in-place...")
    restore_chain_ids(args.pdb_dir, args.out_dir)

if __name__ == "__main__":
    main()

