from Bio import SeqIO
import argparse
import os
import sys
sys.path.insert(0, os.path.abspath("wrapper"))
from wrapper.utils.struct_manager import load_pose_from_pdb

def get_original_sequence(pdb_path):
    pose = load_pose_from_pdb(pdb_path)
    return pose.seq  # str: full sequence, concatenated from all chains

def get_designed_sequence(fasta_path):
    return str(next(SeqIO.parse(fasta_path, "fasta")).seq)

def compare_sequences(orig_seq, designed_seq):
    if len(orig_seq) != len(designed_seq):
        print(f"Sequence length mismatch! Original: {len(orig_seq)}, Designed: {len(designed_seq)}")
        return

    print(f"Comparing sequences of length {len(orig_seq)}...")
    changed = 0
    for i, (o, d) in enumerate(zip(orig_seq, designed_seq)):
        if o != d:
            print(f"Residue {i:3}: {o} to {d}")
            changed += 1

    print(f"\nSummary: {changed} residues changed.")
    if changed == 0:
        print("All positions are preserved no CDR design detected!")
    else:
        print("Framework residues likely preserved; check against your CDR mask.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare original vs. designed FAMPnn sequence.")
    parser.add_argument("--pdb", required=True, help="Path to input PDB file (original)")
    parser.add_argument("--fasta", required=True, help="Path to designed FASTA file from FAMPnn")
    args = parser.parse_args()

    orig = get_original_sequence(args.pdb)
    designed = get_designed_sequence(args.fasta)

    compare_sequences(orig, designed)
