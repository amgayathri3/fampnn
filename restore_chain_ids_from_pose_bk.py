from rfantibody.util.pose import Pose

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def load_fampnn_sequence(pdb_path):
    seq_dict = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                resn = line[17:20].strip()
                resi = int(line[22:26].strip())
                aa = three_to_one.get(resn, 'X')
                if chain not in seq_dict:
                    seq_dict[chain] = []
                if len(seq_dict[chain]) == 0 or resi != seq_dict[chain][-1][0]:
                    seq_dict[chain].append((resi, aa))

    return {
        chain: ''.join(aa for (_, aa) in sorted(entries))
        for chain, entries in seq_dict.items()
    }

def match_chains_by_sequence(original_seq_dict, designed_seq_dict):
    used = set()
    mapping = {}
    for orig_chain, orig_seq in original_seq_dict.items():
        for design_chain, design_seq in designed_seq_dict.items():
            if design_chain in used:
                continue
            if orig_seq[:10] == design_seq[:10]:  # simple heuristic match
                mapping[design_chain] = orig_chain
                used.add(design_chain)
                break
    return mapping

def restore_chains(original_pdb, fampnn_output_pdb, output_pdb):
    original_pose = Pose.from_pdb(original_pdb)
    designed_seqs = load_fampnn_sequence(fampnn_output_pdb)

    mapping = match_chains_by_sequence(original_pose.seq_dict, designed_seqs)

    new_seq_dict = {}
    for design_chain, orig_chain in mapping.items():
        new_seq_dict[orig_chain] = designed_seqs[design_chain]

    output_pose = original_pose.clone()
    output_pose.apply_sequence_changes(new_seq_dict)
    output_pose.to_pdb(output_pdb)
    print(f"Saved with restored chain IDs to {output_pdb}")

