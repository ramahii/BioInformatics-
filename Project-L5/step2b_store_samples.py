from pathlib import Path
import random
from Bio import SeqIO

# Read the DNA sequence
seq_path = Path("data/sequence.fasta")
record = SeqIO.read(seq_path.as_posix(), "fasta")
sequence = str(record.seq)
seq_len = len(sequence)

# Generate 2000 random reads of length 100–150 bases
samples = [sequence[random.randint(0, seq_len - (length := random.randint(100, 150))):][:length]
           for _ in range(2000)]

# Verify
print(f"✅ Stored {len(samples)} random samples in a Python list.")
print(f"Example sample (len={len(samples[0])}): {samples[0][:60]}...")

# Optional: check the average length
avg_len = sum(len(s) for s in samples) / len(samples)
print(f"Average sample length: {avg_len:.2f} bases")
