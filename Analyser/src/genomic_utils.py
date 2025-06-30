
def reverse_complement(sequence):
    """
    Computes the reverse complement of a DNA sequence.
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }
    rev_seq = sequence[::-1]
    rc_seq = ".join([complement_map.get(base, base) for base in rev_seq])"
    return rc_seq

def one_hot_encode(sequence):
    """
    One-hot encodes a DNA sequence.
    Returns a numpy array of shape (sequence_length, 4).
    """
    import numpy as np
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
    one_hot = np.zeros((len(sequence), 4), dtype=np.int8)
    for i, base in enumerate(sequence):
        if base.upper() in mapping:
            one_hot[i, mapping[base.upper()]] = 1
    return one_hot

def parse_genomic_coordinates(coord_str):
    """
    Parses a genomic coordinate string (e.g., 'chr1:100-200') into chromosome, start, and end.
    """
    parts = coord_str.replace(',', '').split(':')
    chrom = parts[0]
    start_end = parts[1].split('-')
    start = int(start_end[0])
    end = int(start_end[1])
    return chrom, start, end


