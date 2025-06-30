import os
import urllib.request

def download_file(url, dest_folder):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    file_name = os.path.join(dest_folder, os.path.basename(url))
    print(f"Downloading {url} to {file_name}...")
    urllib.request.urlretrieve(url, file_name)
    print("Download complete.")

if __name__ == "__main__":
    chm13_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
    dest_folder = "../data/chm13_reference"
    download_file(chm13_url, dest_folder)


