# Human_Genome_Loader 🚀

A lightweight, efficient toolkit for loading and processing genomic data (FASTA, 2bit, BigWig, BedGraph, etc.) tailored for deep learning workflows.

---

## Table of Contents

1. [Overview](#overview)  
2. [Features](#features)  
3. [Installation](#installation)  
4. [Usage Examples](#usage-examples)  
   - [Loading DNA sequences](#loading-dna-sequences)  
   - [Processing signal tracks (BigWig/BedGraph)](#processing-signal-tracks)  
   - [Using with Keras/PyTorch generators](#using-with-keraspytorch-generators)  
5. [Dependencies](#dependencies)  
6. [Development & Contribution](#development--contribution)  
7. [Licensing](#licensing)  
8. [Citation](#citation)

---

## Overview

**Human_Genome_Loader** is designed for deep-learning projects involving genomic data. It offers:

- Support for popular formats: FASTA, 2bit, BigWig, BedGraph, BED.
- On-the-fly, memory-efficient data loading.
- Pythonic APIs and dataset generators compatible with TensorFlow (Keras) and PyTorch.
- Minimal external dependencies.
- Reproducible training via configurable random seeds.

---

## Features

- **Efficient data streaming** from disk  
- **2bit support** for reduced storage and faster access  
- **BigWig / BedGraph wrappers** for signal extraction  
- **BED region selection** with interval utilities  
- **Generators** to feed sequenced or signal-labeled batches into DL models  
- **Flexible and extensible**—easily add new formats  

---

## Installation

```bash
git clone https://github.com/Ftyigffifygf/Human_Genome_Loader.git
cd Human_Genome_Loader
python setup.py develop
```

> ⚠️ Using `develop` is ideal during active development to reflect code changes immediately.

---

## Usage Examples

### Loading DNA sequences (2bit/FASTA)

```python
from hgl.wrapper import TwoBitWrapper

twobit = TwoBitWrapper("hg38.2bit")
sequence = twobit.fetch(seq_name="chr12", start=100000, end=100100)
print(sequence)  # prints 100 bp DNA sequence
```

### Processing signal tracks (BigWig/BedGraph)

```python
from hgl.wrapper import BigWigWrapper

bw = BigWigWrapper("sample_signal.bw")
signal = bw.fetch(seq_name="chr1", start=20000, end=20100)
print(signal)  # a numpy array of signal values
```

### Using with Keras/PyTorch generators

```python
from keras.models import Sequential
from keras.layers import Conv1D, Flatten, Dense
from hgl.generator import BedGraphGenerator

# Set up wrappers
twobit = TwoBitWrapper("hg38.2bit")
bg = BedGraphWrapper("intervals.bedGraph.gz")

# Create generator
gen = BedGraphGenerator(bg, twobit, seq_len=101, batch_size=64, shuffle=True)

# Build Keras model
model = Sequential([
    Conv1D(32, 21, input_shape=(101, 4)),
    Flatten(),
    Dense(1, activation="sigmoid")
])
model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])

# Train
model.fit(gen, steps_per_epoch=200, epochs=10)
```

_PyTorch Dataset/Loader support coming soon!_

---

## Dependencies

### Required

- `pyfaidx` – FASTA/2bit sequence indexing  
- `py2bit` – 2bit file support  
- `pyBigWig` – BigWig/BigBed support  
- `pybedtools` – BED file utilities  
- `quicksect` – Interval tree querying  
- `numpy`, `pandas`, `scikit-learn`, `tqdm` – data handling & progress bars  

### Optional

- `biopython` – for compressed FASTA access  
- `tensorflow/keras` or `torch` – model training frameworks  

> Using Anaconda (Python 3.8+) is recommended for smoother dependency handling.

---

## Development & Contribution

- Clone the repo and install with `python setup.py develop`
- Add new data format support or generator backends
- Keep compatibility up-to-date
- Ensure low memory usage (aim <2GB)
- Write tests for new modules
- Follow existing code style and seed initializers for reproducibility

---

## Licencing

Distributed under the **MIT License**. See `LICENSE` for details.

---

## Citation

If utilizing this tool in publications, please cite:

```
COMING SOON — citation details to be added
```

---


---

### Next Steps

- [ ] Add PyTorch `Dataset` / `DataLoader` support  
- [ ] Expand format support (e.g., BAM, CRAM)  
- [ ] Add multi-task generator support  
- [ ] Improve documentation with notebooks & tutorials
