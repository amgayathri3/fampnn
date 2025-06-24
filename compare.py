from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
from Bio.Seq import Seq

# Load input PDB and extract sequence
p = PDBParser(QUIET=True)
structure = p.get_structure("input", "6a3w_rbd_G496T500G502_knob_3_1_0.pdb")

seqs = []
for model in structure:
    for chain in model:
        ppb = PPBuilder()
        for pp in ppb.build_peptides(chain):
            seqs.append(str(pp.get_sequence()))

input_seq = "".join(seqs)

# Load output FASTA sequence
out_fasta = "outputs/cdr_test_run/fastas/6a3w_rbd_G496T500G502_knob_3_1_0_sample0.fasta"
out_seq = str(list(SeqIO.parse(out_fasta, "fasta"))[0].seq)

# Compare and report mutations
for i, (a, b) in enumerate(zip(input_seq, out_seq)):
    if a != b:
        print("Mutation at position {}: {} to {}".format(i + 1, a, b))
