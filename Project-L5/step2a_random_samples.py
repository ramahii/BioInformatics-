from pathlib import Path
import random
from Bio import SeqIO

# Load the original sequence
SEQ_PATH = Path("data/sequence.fasta")
record = SeqIO.read(SEQ_PATH.as_posix(), "fasta")
sequence = str(record.seq)
seq_len = len(sequence)

# Generate 2000 random samples (reads) of length 100–150
samples = []
for _ in range(2000):
    length = random.randint(100, 150)
    start = random.randint(0, seq_len - length)
    sample = sequence[start:start + length]
    samples.append(sample)

print(f"✅ Generated {len(samples)} samples")
print("Example sample:", samples[0][:60] + "...")

# (optional) save to file so we can inspect later
out_path = Path("data/random_samples.txt")
with out_path.open("w") as f:
    for s in samples:
        f.write(s + "\n")

print(f"Samples saved to {out_path}")
