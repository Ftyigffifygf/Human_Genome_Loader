from Bio.Seq import Seq

# Bonus: DNA → RNA → Protein simulation
sample_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print("DNA :", sample_dna)

rna = sample_dna.transcribe()
print("RNA :", rna)

protein = rna.translate()
print("Protein:", protein)
