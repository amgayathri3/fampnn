import os
import torch

from rfantibody.rf2.modules.pose_util import pose_from_remarked, Pose, pose_to_remarked_pdblines, pdblines_to_pdb
from rfantibody.rf2.modules.util import parse_pdb_coords_seq

def restore_chains(original_pdb_path, modified_pdb_path, output_pdb_path):
    """
    Restore HLT chain IDs and REMARK annotations to a FAMPnn-generated PDB.
    
    Args:
        original_pdb_path (str): Path to the original HLT input PDB with REMARK annotations.
        modified_pdb_path (str): Path to the FAMPnn output PDB (no HLT labels).
        output_pdb_path (str): Path to write the restored HLT-labeled output PDB.
    """
    # Load original pose with full metadata
    pose = pose_from_remarked(original_pdb_path)

    # Load new coordinates and sequence from modified PDB
    new_xyz, new_seq = parse_pdb_coords_seq(modified_pdb_path)

    # Build new pose object using original metadata and new output
    restored_pose = Pose(
        xyz=new_xyz,
        seq=new_seq,
        atom_mask=pose.atom_mask,
        cdrs=pose.cdrs,
        idx=pose.idx,
        chain_dict=pose.chain_dict
    )

    # Re-generate PDB lines with HLT chain IDs and REMARKs
    pdblines = pose_to_remarked_pdblines(restored_pose)

    # Write final restored PDB
    pdblines_to_pdb(pdblines, output_pdb_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pose", required=True, help="Path to original REMARKED HLT PDB")
    parser.add_argument("--samples_dir", required=True, help="Directory with FAMPnn sample PDBs")
    parser.add_argument("--tag", required=True, help="Prefix (e.g., 6a3w_rbd) to match sample PDB files")

    args = parser.parse_args()

    # Load once
    pose = pose_from_remarked(args.input_pose)

    for fname in os.listdir(args.samples_dir):
        if fname.startswith(args.tag) and fname.endswith(".pdb"):
            sample_path = os.path.join(args.samples_dir, fname)
            print(f"[INFO] Restoring {fname}...")

            try:
                new_xyz, new_seq = parse_pdb_coords_seq(sample_path)
                new_pose = Pose(
                    xyz=new_xyz,
                    seq=new_seq,
                    atom_mask=pose.atom_mask,
                    cdrs=pose.cdrs,
                    idx=pose.idx,
                    chain_dict=pose.chain_dict
                )
                pdblines = pose_to_remarked_pdblines(new_pose)
                pdblines_to_pdb(pdblines, sample_path)
            except Exception as e:
                print(f"[ERROR] Failed to restore {fname}: {e}")

    print("[DONE] All files restored with original chain IDs.")

