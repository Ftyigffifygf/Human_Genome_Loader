from genomekit.fasta_handler import FastaHandler
from genomekit.reference_indexer import ensure_index
from genomekit.variant_loader import VariantLoader
from genomekit.annotator import GenomeAnnotator
from genomekit.crispr_editor import CRISPREditor

import sys

def run(fasta_path, variant_path, chrom, start, end, gRNA):
    ensure_index(fasta_path)
    fasta = FastaHandler(fasta_path)
    variant = VariantLoader(variant_path)
    annotator = GenomeAnnotator(fasta, variant)
    editor = CRISPREditor()

    seq, vars = annotator.annotate_region(chrom, int(start), int(end))
    print(f"Original Sequence: {seq}")
    print(f"Variants: {vars}")

    edits = editor.simulate_cut(seq, gRNA)
    print("\n--- Simulated Edits ---")
    for pos, edit_seq in edits:
        print(f"Cut at: {pos}, Edited Sequence: {edit_seq}")

if __name__ == "__main__":
    run(*sys.argv[1:])
