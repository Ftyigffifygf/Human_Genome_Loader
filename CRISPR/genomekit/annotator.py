class GenomeAnnotator:
    def __init__(self, fasta_handler, variant_loader):
        self.fasta_handler = fasta_handler
        self.variant_loader = variant_loader

    def annotate_region(self, chrom, start, end):
        seq = self.fasta_handler.get_sequence(chrom, start, end)
        variants = self.variant_loader.get_variants(chrom, start, end)
        return seq, variants
