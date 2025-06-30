
from src.data_loader import GenomeDataLoader
from src.genomic_utils import parse_genomic_coordinates
import os

# Path to the downloaded CHM13 reference genome
# Make sure to run scripts/download_chm13.py first
chm13_ref_path = os.path.join(os.path.dirname(__file__), "..", "data", "chm13_reference", "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")

if not os.path.exists(chm13_ref_path):
    print(f"Error: CHM13 reference genome not found at {chm13_ref_path}")
    print("Please run `python scripts/download_chm13.py` first.")
else:
    loader = GenomeDataLoader(chm13_ref_path)

    # Example 1: Load a sequence from a specific region
    coord_str = "chr1:10000-10100"
    chrom, start, end = parse_genomic_coordinates(coord_str)
    sequence = loader.get_sequence(chrom, start, end)
    print(f"\nSequence from {coord_str}:\n{sequence}")

    # Example 2: Load a sequence and reverse complement it
    from src.genomic_utils import reverse_complement
    rc_sequence = reverse_complement(sequence)
    print(f"\nReverse complement of the sequence:\n{rc_sequence}")

    # Example 3: One-hot encode the sequence
    from src.genomic_utils import one_hot_encode
    one_hot_seq = one_hot_encode(sequence)
    print(f"\nOne-hot encoded sequence (first 5 bases):\n{one_hot_seq[:5]}")


