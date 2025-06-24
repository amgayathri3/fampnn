import sys
import os
sys.path.insert(0, os.path.abspath("wrapper"))

from wrapper.utils.struct_manager import load_pose_from_pdb

pose = load_pose_from_pdb("6a3w_rbd_G496T500G502_knob_3_1_0.pdb")
print("CDR Dict:", pose.cdr_dict)
print("Chain Lengths:", pose.chain_length_dict)

