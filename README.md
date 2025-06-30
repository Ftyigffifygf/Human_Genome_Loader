# 🧬 Human Genome Loader (Python Project)

This project downloads and loads the **entire human genome**, chromosome-by-chromosome, from NCBI using Python and Biopython.

## 🔧 Requirements

Install required libraries:

```bash
pip install biopython tqdm matplotlib
```

## 📁 Files

- `main.py`: Main script to download, parse, and analyze genome data.
- `genome_data/`: Directory where `.fasta` chromosome files are saved.

## 🚀 How to Run

```bash
python main.py
```

It will:
1. Download all 23 pairs human chromosomes (from NCBI).
2. Save them as FASTA files.
3. Load them into memory.
4. Print total genome size.
5. Plot GC content across chromosomes.

## 📈 Example Output

- Chromosome GC content bar plot
- Total genome length printed

## 📬 Note

Set your email in `Entrez.email = "your_email@example.com"` so NCBI knows who is accessing.

---

Inspired by OpenAI | Powered by Biopython 🧪
