
import pyfaidx
import pysam

class GenomeDataLoader:
    """
    A data loader for handling genomic data, including reference sequences and variants.
    """

    def __init__(self, reference_path):
        """
        Initializes the data loader with a path to a reference genome file (FASTA or 2bit).
        `pyfaidx` is used for efficient sequence retrieval.
        """
        self.reference = pyfaidx.Fasta(reference_path)

    def get_sequence(self, chrom, start, end):
        """
        Extracts a genomic sequence from the reference genome.
        Coordinates are 1-based.
        """
        return self.reference[chrom][start-1:end].seq

    def apply_variants(self, chrom, start, end, vcf_file):
        """
        Applies variants from a VCF file to a reference sequence to generate a personalized sequence.
        This is a simplified implementation for demonstration purposes.
        """
        sequence = list(self.get_sequence(chrom, start, end))
        vcf = pysam.VariantFile(vcf_file)

        for record in vcf.fetch(chrom, start-1, end):
            pos_in_sequence = record.pos - start
            if 0 <= pos_in_sequence < len(sequence):
                # This example only handles simple SNPs
                if len(record.ref) == 1 and len(record.alts[0]) == 1:
                    sequence[pos_in_sequence] = record.alts[0]

        return "".join(sequence)


