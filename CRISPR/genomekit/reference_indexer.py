import os
import pysam

def ensure_index(fasta_path):
    if not os.path.exists(fasta_path + ".fai"):
        print("Indexing reference FASTA...")
        pysam.faidx(fasta_path)
