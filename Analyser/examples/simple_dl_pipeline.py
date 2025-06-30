import numpy as np
from src.data_loader import GenomeDataLoader
from src.genomic_utils import one_hot_encode, parse_genomic_coordinates
import os

# This is a conceptual example. A real deep learning setup would involve
# a framework like TensorFlow or PyTorch and more complex data handling.

# Path to the downloaded CHM13 reference genome
chm13_ref_path = os.path.join(os.path.dirname(__file__), "..", "data", "chm13_reference", "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")

if not os.path.exists(chm13_ref_path):
    print(f"Error: CHM13 reference genome not found at {chm13_ref_path}")
    print("Please run `python scripts/download_chm13.py` first.")
else:
    loader = GenomeDataLoader(chm13_ref_path)

    # Define a genomic region to simulate data loading for a DL model
    coord_str = "chr1:10000-10200"
    chrom, start, end = parse_genomic_coordinates(coord_str)

    # Simulate loading a batch of genomic sequences
    batch_size = 4
    sequence_length = end - start
    
    # In a real scenario, you would load different sequences for each item in the batch
    # Here, we'll just use the same sequence for simplicity.
    base_sequence = loader.get_sequence(chrom, start, end)
    
    # One-hot encode the base sequence
    one_hot_base_sequence = one_hot_encode(base_sequence)

    # Create a batch of one-hot encoded sequences
    # Shape: (batch_size, sequence_length, 4)
    data_batch = np.array([one_hot_base_sequence] * batch_size)

    print(f"\nSimulated data batch shape for DL model: {data_batch.shape}")
    print(f"Example of one-hot encoded sequence (first 5 bases of first item in batch):\n{data_batch[0, :5]}")

    # Conceptual placeholder for a deep learning model
    # In a real application, you would define and train your model here.
    def conceptual_dl_model(input_data):
        print("\nConceptual DL model received data.")
        print(f"Input data shape: {input_data.shape}")
        # Simulate some processing
        output = np.sum(input_data, axis=(1, 2))
        return output

    # Pass the data batch to the conceptual model
    model_output = conceptual_dl_model(data_batch)
    print(f"Conceptual DL model output (sum of bases per sequence): {model_output}")

    print("\nThis example demonstrates how the data_loader and genomic_utils can prepare data for a deep learning pipeline.")


