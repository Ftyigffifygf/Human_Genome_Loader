from Bio.Seq import Seq
import re

class CRISPREditor:
    def __init__(self, pam_seq="NGG"):
        self.pam_pattern = pam_seq.replace("N", "[ATCG]")

    def find_targets(self, sequence, guide_rna):
        pattern = guide_rna + self.pam_pattern
        matches = [(m.start(), m.end()) for m in re.finditer(pattern, sequence)]
        return matches

    def simulate_cut(self, sequence, guide_rna):
        targets = self.find_targets(sequence, guide_rna)
        edits = []
        for start, end in targets:
            cut_pos = start + len(guide_rna) - 3
            edited = sequence[:cut_pos] + '-' + sequence[cut_pos+1:]
            edits.append((cut_pos, edited))
        return edits
