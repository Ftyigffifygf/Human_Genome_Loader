
from src.data_loader import GenomeDataLoader
from src.genomic_utils import parse_genomic_coordinates
import os

# Path to the downloaded CHM13 reference genome
# Make sure to run scripts/download_chm13.py first
chm13_ref_path = os.path.join(os.path.dirname(__file__), "..", "data", "chm13_reference", "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")

# Create a dummy VCF file for demonstration
dummy_vcf_path = os.path.join(os.path.dirname(__file__), "dummy.vcf")
with open(dummy_vcf_path, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##contig=<ID=chr1,length=248956422>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    f.write("chr1\t10005\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n")
    f.write("chr1\t10050\t.\tG\tC\t.\tPASS\t.\tGT\t1/1\n")

if not os.path.exists(chm13_ref_path):
    print(f"Error: CHM13 reference genome not found at {chm13_ref_path}")
    print("Please run `python scripts/download_chm13.py` first.")
else:
    loader = GenomeDataLoader(chm13_ref_path)

    # Example: Apply variants to a sequence
    coord_str = "chr1:10000-10100"
    chrom, start, end = parse_genomic_coordinates(coord_str)

    original_sequence = loader.get_sequence(chrom, start, end)
    print(f"\nOriginal sequence from {coord_str}:\n{original_sequence}")

    personalized_sequence = loader.apply_variants(chrom, start, end, dummy_vcf_path)
    print(f"\nPersonalized sequence with variants:\n{personalized_sequence}")

# Clean up dummy VCF file
os.remove(dummy_vcf_path)


