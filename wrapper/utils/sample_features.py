from rfantibody.util.pose import Pose

class SampleFeatures:
    def __init__(self, pose: Pose):
        self.pose = pose
        self.fixed_dict = {}

    def set_fixed_from_cdrs(self):
        for chain, seq in self.pose.seq_dict.items():
            fixed = [True] * len(seq)
            for start, end in self.pose.cdr_dict.get(chain, []):
                for i in range(start, end + 1):
                    fixed[i] = False
            self.fixed_dict[chain] = fixed

    def to_model_input(self):
        return {
            "coords": self.pose.coords_array(),
            "mask": self.pose.mask_array(),
            "seq_dict": self.pose.seq_dict,
            "chains": list(self.pose.seq_dict.keys()),
        }

