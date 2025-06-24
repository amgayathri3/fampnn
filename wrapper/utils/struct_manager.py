from rfantibody.util.pose import Pose

def load_pose_from_pdb(pdb_path):
    return Pose.from_pdb(pdb_path)

