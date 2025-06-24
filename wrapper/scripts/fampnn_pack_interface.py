import os
import argparse
import glob
import sys
from fampnn.inference.pack import main as fampnn_pack_main

def run_pack_wrapper(
    pdb_path=None,
    pdb_dir=None,
    out_dir=None,
    checkpoint_path=None,
    fixed_pos_csv=None
):
    # Setup input PDBs
    input_dir = os.path.join(out_dir, "input_pdbs")
    os.makedirs(input_dir, exist_ok=True)

    pdb_files = []
    if pdb_dir:
        pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    elif pdb_path:
        pdb_files = [pdb_path]
    else:
        raise ValueError("Provide either --pdb_path or --pdb_dir.")

    if not pdb_files:
        raise FileNotFoundError("No PDB files found.")

    # Copy pdbs and write key list (keep full filenames)
    keylist_path = os.path.join(out_dir, "temp_keys.txt")
    with open(keylist_path, "w") as f:
        for pdb in pdb_files:
            basename = os.path.basename(pdb)
            target_path = os.path.join(input_dir, basename)
            os.system(f"cp {pdb} {target_path}")
            f.write(basename + "\n")

    # Simulate CLI args for pack.py
    sys.argv = [
        "pack.py",
        f"checkpoint_path={checkpoint_path}",
        f"pdb_dir={input_dir}",
        f"pdb_key_list={keylist_path}",
        f"out_dir={out_dir}"
    ]
    if fixed_pos_csv:
        sys.argv.append(f"fixed_pos_csv={fixed_pos_csv}")

    fampnn_pack_main()

def main():
    parser = argparse.ArgumentParser(description="FAMPnn Wrapper for Sidechain Packing")
    parser.add_argument("--pdb_path", help="Path to a single PDB file")
    parser.add_argument("--pdb_dir", help="Path to a folder of PDBs")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--checkpoint_path", required=True, help="Path to model weights")
    parser.add_argument("--fixed_pos_csv", help="Optional fixed positions CSV")

    args = parser.parse_args()

    run_pack_wrapper(
        pdb_path=args.pdb_path,
        pdb_dir=args.pdb_dir,
        out_dir=args.out_dir,
        checkpoint_path=args.checkpoint_path,
        fixed_pos_csv=args.fixed_pos_csv,
    )

if __name__ == "__main__":
    main()

