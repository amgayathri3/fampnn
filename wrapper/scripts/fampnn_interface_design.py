import os
import sys
import argparse
import torch
import time
import json

from wrapper.utils.sample_features import SampleFeatures
from wrapper.utils.io import write_fasta
from wrapper.utils.pose import Pose
from inference.seq_design import run_seq_design

def main():
    parser = argparse.ArgumentParser(description="Run FAMPnn design with RFantibody-style masking.")
    parser.add_argument("--pdb_path", type=str, required=True, help="Path to input PDB file")
    parser.add_argument("--checkpoint_path", type=str, required=True, help="Path to FAMPnn model weights (.pt)")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save output files")
    parser.add_argument("--num_seqs", type=int, default=8, help="Number of sequences to generate")
    parser.add_argument("--device", type=str, default="cuda", help="Device to run on")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pdb_name = os.path.splitext(os.path.basename(args.pdb_path))[0]

    pose = Pose(args.pdb_path)

    # Construct fixed_positions_dict using CDR logic
    fixed_positions_dict = {}
    for chain, seq in pose.seq_dict.items():
        fixed_positions_dict[chain] = []
        for idx, aa in enumerate(seq):
            abs_res = pose.get_absolute_residue_index(chain, idx)
            if not pose.cdr_dict.get(chain, {}).get(abs_res, False):
                fixed_positions_dict[chain].append(abs_res)

    # Prepare features
    sample = SampleFeatures(
        pdb_path=args.pdb_path,
        fixed_positions_dict=fixed_positions_dict,
        cdr_dict=pose.cdr_dict,
    )

    # Run design
    designed_seqs = run_seq_design(
        features=sample,
        checkpoint_path=args.checkpoint_path,
        num_seqs=args.num_seqs,
        device=args.device
    )

    # Write output
    fasta_out = os.path.join(args.out_dir, f"{pdb_name}_fampnn_design.fasta")
    write_fasta(designed_seqs, fasta_out)

    print(f"Designs saved to: {fasta_out}")

if __name__ == "__main__":
    main()

