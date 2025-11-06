from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import os
import requests

# Step 1: Influenza genomes with NCBI accessions
genomes = {
    "Influenza_A_H1N1": "NC_002023.1",
    "Influenza_A_H3N2": "NC_007366.1",
    "Influenza_A_H5N1": "NC_007366.2",
    "Influenza_A_H7N9": "NC_026433.1",
    "Influenza_B_Yamagata": "NC_002205.1",
    "Influenza_B_Victoria": "NC_002204.1",
    "Influenza_C": "NC_006308.1",
    "Influenza_A_H9N2": "NC_011516.1",
    "Influenza_A_H10N8": "NC_026432.1",
    "Influenza_A_H1N2": "NC_007373.1"
}

os.makedirs("genomes", exist_ok=True)

# Step 2: Download proper FASTA using NCBI efetch
for name, accession in genomes.items():
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    response = requests.get(url)
    with open(f"genomes/{name}.fasta", "w", encoding="utf-8") as f:
        f.write(response.text)

# Step 3: Define EcoRI restriction site
ECO_RI_SITE = "GAATTC"

# Step 4: Parse genomes and apply restriction digestion
fragment_data = {}

for name in genomes.keys():
    fasta_path = f"genomes/{name}.fasta"
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            print(f"[Warning] No valid sequences found in {name}. Skipping.")
            continue
        sequence = str(records[0].seq)
        fragments = sequence.split(ECO_RI_SITE)
        fragment_lengths = [len(frag) for frag in fragments if len(frag) > 0]
        fragment_data[name] = fragment_lengths
    except Exception as e:
        print(f"[Error] Could not parse {name}: {e}")

# Step 5: Find genome with most EcoRI fragments
if not fragment_data:
    raise ValueError("No valid genomes could be processed. Check internet connection or accessions.")

most_fragments_name = max(fragment_data, key=lambda x: len(fragment_data[x]))
print(f"Genome with the most EcoRI fragments: {most_fragments_name} ({len(fragment_data[most_fragments_name])} fragments)")

# Step 6: Function to simulate and plot gel
def plot_gel(fragment_lengths, title, filename):
    min_bp, max_bp = min(fragment_lengths), max(fragment_lengths)
    positions = [np.interp(np.log(size + 1), [np.log(min_bp + 1), np.log(max_bp + 1)], [0.9, 0.1]) for size in fragment_lengths]
    sorted_pairs = sorted(zip(fragment_lengths, positions), key=lambda x: x[0], reverse=True)
    sorted_sizes, sorted_positions = zip(*sorted_pairs)

    fig, ax = plt.subplots(figsize=(3, 6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    ax.add_patch(plt.Rectangle((0.1, 0), 0.8, 1, color='black'))
    ax.plot([0.1, 0.9, 0.9, 0.1, 0.1], [0, 0, 1, 1, 0], color='white', linewidth=2)

    for size, pos in zip(sorted_sizes, sorted_positions):
        ax.hlines(pos, 0.2, 0.8, color='white', linewidth=3)

    ax.text(0.02, 0.9, "3000 bp -", color='black', fontsize=9)
    ax.text(0.02, 0.5, "1500 bp -", color='black', fontsize=9)
    ax.text(0.02, 0.2, "500 bp -", color='black', fontsize=9)
    plt.title(title, color='black', pad=20)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# Step 7: Generate and save gels
os.makedirs("gels", exist_ok=True)

for name, lengths in fragment_data.items():
    plot_gel(lengths, f"{name} - EcoRI digest", f"gels/{name}.png")

print("Gel simulations saved in the 'gels' folder.")
