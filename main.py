from Bio import Entrez, SeqIO
from tqdm import tqdm
import os
import matplotlib.pyplot as plt

Entrez.email = "your_email@example.com"

chromosomes = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10"
}

def download_chromosomes(chrom_list, out_dir="genome_data"):
    os.makedirs(out_dir, exist_ok=True)
    for name, ncbi_id in tqdm(chrom_list.items(), desc="Downloading chromosomes"):
        filename = f"{out_dir}/chr{name}.fasta"
        if os.path.exists(filename):
            print(f"✔️ Chromosome {name} already downloaded.")
            continue
        try:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
            with open(filename, "w") as out_f:
                out_f.write(handle.read())
            handle.close()
        except Exception as e:
            print(f"❌ Error downloading chr{name}: {e}")

def load_all_sequences(folder="genome_data"):
    genome = {}
    for file in os.listdir(folder):
        if file.endswith(".fasta"):
            path = os.path.join(folder, file)
            record = SeqIO.read(path, "fasta")
            chrom = file.replace("chr", "").replace(".fasta", "")
            genome[chrom] = record.seq
    return genome

def plot_gc_content(genome):
    gc_contents = {}
    for chrom, seq in genome.items():
        g = seq.count("G")
        c = seq.count("C")
        gc_contents[chrom] = 100 * (g + c) / len(seq)

    plt.bar(gc_contents.keys(), gc_contents.values())
    plt.title("GC Content per Chromosome")
    plt.ylabel("GC%")
    plt.xlabel("Chromosome")
    plt.show()

if __name__ == "__main__":
    download_chromosomes(chromosomes)
    genome = load_all_sequences()
    print(f"Total genome size: {sum(len(seq) for seq in genome.values()):,} base pairs")
    plot_gc_content(genome)
