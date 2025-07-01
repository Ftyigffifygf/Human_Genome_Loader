# Human_Genome_Loader 🚀

A lightweight, efficient toolkit for loading and processing genomic data (FASTA, 2bit, BigWig, BedGraph, BED, etc.) tailored for deep-learning and data-science workflows.

---

## Table of Contents

- [About](#about)  
- [Features](#features)  
- [Installation](#installation)  
- [Quick Start](#quick-start)  
- [API Reference](#api-reference)  
  - [Sequence Wrappers](#sequence-wrappers)  
  - [Signal Wrappers](#signal-wrappers)  
  - [Generators](#generators)  
- [Command-Line Utilities](#command-line-utilities)  
- [Development & Testing](#development-testing)  
- [Contributing](#contributing)  
- [License](#license)  

---

## About

Human_Genome_Loader provides:

- On-the-fly, memory-efficient streaming of genomic files.  
- Pythonic interfaces for popular formats: FASTA, 2bit, BigWig, BedGraph, BED.  
- Dataset generators compatible with TensorFlow/Keras and PyTorch (coming soon).  
- Minimal dependencies and reproducible behavior via configurable seeds.  

---

## Features

- Efficient random access of large genome assemblies (2bit/FASTA)  
- BigWig & BedGraph wrappers for signal extraction as NumPy arrays  
- BED‐based interval selection and fast interval queries  
- Keras `Sequence` & custom PyTorch `Dataset` generators  
- Shell installation script and example pipelines  
- Extensible design—add new formats or generators  

---

## Installation

Clone and install in developer mode (reflects code changes immediately):

```bash
git clone https://github.com/Ftyigffifygf/Human_Genome_Loader.git
cd Human_Genome_Loader
bash install.sh
```

Alternatively, install manually:

```bash
pip install pyfaidx py2bit pyBigWig pybedtools quicksect numpy pandas scikit-learn tqdm
pip install .  # from project root
```

> Recommended: Use Python 3.8+ in an Anaconda environment for smoother dependency handling.

---

## Quick Start

### 1. Load DNA Sequences (2bit/FASTA)

```python
from hgl.wrapper import TwoBitWrapper, FastaWrapper

# 2bit example
tb = TwoBitWrapper("hg38.2bit")
seq = tb.fetch(seq_name="chr12", start=100_000, end=100_100)
print(seq)  # 100-bp DNA string

# FASTA example
fa = FastaWrapper("chr1.fa.gz")
seq2 = fa.fetch("chr1", 200_000, 200_050)
print(seq2)
```

### 2. Extract Signal Tracks (BigWig/BedGraph)

```python
from hgl.wrapper import BigWigWrapper, BedGraphWrapper

bw = BigWigWrapper("signal.bw")
signal_bw = bw.fetch("chr1", 1_000, 1_100)

bg = BedGraphWrapper("signal.bedGraph.gz")
signal_bg = bg.fetch("chr2", 50_000, 50_100)

print(signal_bw.shape, signal_bg.shape)  # (100,), (100,)
```

### 3. Train with Keras Generator

```python
from hgl.wrapper import TwoBitWrapper, BedGraphWrapper
from hgl.generator import BedGraphGenerator
from keras.models import Sequential
from keras.layers import Conv1D, Flatten, Dense

tb = TwoBitWrapper("hg38.2bit")
bg = BedGraphWrapper("peaks.bedGraph.gz")

gen = BedGraphGenerator(
    bg_wrapper=bg, seq_wrapper=tb,
    seq_len=101, batch_size=32, shuffle=True
)

model = Sequential([
    Conv1D(32, 21, input_shape=(101, 4)),
    Flatten(),
    Dense(1, activation="sigmoid")
])

model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])
model.fit(gen, steps_per_epoch=100, epochs=5)
```

---

## API Reference

### Sequence Wrappers

| Class             | File               | Description                                      |
|-------------------|--------------------|--------------------------------------------------|
| TwoBitWrapper     | `wrapper.py`       | Random access to `.2bit` genome assemblies       |
| FastaWrapper      | `wrapper.py`       | Indexed FASTA (supports gzip/compression)        |

### Signal Wrappers

| Class             | File               | Description                                      |
|-------------------|--------------------|--------------------------------------------------|
| BigWigWrapper     | `wrapper.py`       | Fetches signal from BigWig/BigBed                 |
| BedGraphWrapper   | `wrapper.py`       | Parses BedGraph files (gzip OK)                  |
| BedWrapper        | `wrapper.py`       | Interval queries on BED files                    |

### Generators

| Class               | File               | Framework      | Purpose                                  |
|---------------------|--------------------|----------------|------------------------------------------|
| BedGraphGenerator   | `generator.py`     | Keras (`Sequence`) | Batch signal + one-hot DNA generator   |
| SequenceGenerator   | `generator.py`     | Keras/PyTorch  | Customizable sequence + label generator  |

---

## Command-Line Utilities

- **`main.py`**: Example pipeline combining wrappers & generators.  
- **`bonus_dna_to_protein.py`**: Translates DNA sequences into proteins using codon tables.  
- **`install.sh`**: Automated setup script.  
- **`genomic_pipeline.py`**: Demonstrates full data preprocessing and saving.

Run any script:

```bash
python main.py --genome hg38.2bit --signal peaks.bw --out results/
```

---

## Development & Testing

1. Create a virtual environment and install in develop mode:  
   ```bash
   python setup.py develop
   ```
2. Add new wrappers or generators in `wrapper.py`/`generator.py`.  
3. Write unit tests under `tests/` and run via:
   ```bash
   pytest --maxfail=1 --disable-warnings -q
   ```
4. Maintain memory usage <2 GB in large benchmarks.  
5. Follow PEP8 style and seed RNGs for reproducibility.

---

## Contributing

We welcome improvements! Please:

1. Fork the repo and create a feature branch.  
2. Write tests for new functionality.  
3. Submit a pull request with a clear description.  
4. Ensure CI passes before merging.

---

## License

Distributed under the MIT License. See [LICENSE.txt](LICENSE.txt) for details.
