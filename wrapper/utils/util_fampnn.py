import os
from fampnn.inference import seq_design


def set_global_fixed_dict(fixed_dict):
    """
    Sets the global in-memory fixed_dict used by seq_design.py
    """
    import fampnn.sampling_utils as sampling_utils
    sampling_utils.global_fixed_dict = fixed_dict


def run_fampnn(
    checkpoint_path,
    pdb_dir,
    pdb_key_list,
    out_dir,
    num_seqs=1,
    temperature=1.0,
):
    """
    Runs the FAMPnn design script with in-memory fixed_dict via Hydra overrides
    """

    overrides = [
        f"checkpoint_path={checkpoint_path}",
        f"pdb_dir={pdb_dir}",
        f"pdb_key_list={pdb_key_list}",
        f"out_dir={out_dir}",
        f"num_seqs_per_pdb={num_seqs}",
        f"batch_size=1",
        f"temperature={temperature}",
        f"seed=42",
    ]

    print("Running FAMPnn via in-process Hydra main():")
    print("Overrides:")
    for line in overrides:
        print(line)

    seq_design.main(overrides=overrides)

