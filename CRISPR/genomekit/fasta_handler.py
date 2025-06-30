import pysam

class FastaHandler:
    def __init__(self, fasta_path):
        self.fasta = pysam.FastaFile(fasta_path)

    def get_sequence(self, chrom, start, end):
        return self.fasta.fetch(chrom, start, end)
