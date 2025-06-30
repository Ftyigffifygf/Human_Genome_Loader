import h5py
import numpy as np

class VariantLoader:
    def __init__(self, hdf5_file):
        self.data = h5py.File(hdf5_file, 'r')

    def get_variants(self, chrom, start, end):
        chrom_data = self.data[chrom]
        start_pos = chrom_data['pos'][:]
        mask = (start_pos >= start) & (start_pos <= end)
        return chrom_data['variants'][mask]
