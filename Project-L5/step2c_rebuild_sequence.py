from pathlib import Path
from Bio import SeqIO

seq_path = Path("data/sequence.fasta")
record = SeqIO.read(seq_path.as_posix(), "fasta")
original_seq = str(record.seq)

samples_path = Path("data/random_samples.txt")
samples = [line.strip() for line in samples_path.open() if line.strip()]

print(f"Loaded {len(samples)} random samples.")

def overlap(a: str, b: str, min_overlap: int = 20) -> int:
   
    start = 0
    while True:
        start = a.find(b[:min_overlap], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def assemble(reads: list[str]) -> str:
    reads = reads.copy()
    assembled = reads.pop(0)
    while reads:
        best_olen = 0
        best_read = None
        best_i = None
        for i, r in enumerate(reads):
            olen = overlap(assembled, r)
            if olen > best_olen:
                best_olen = olen
                best_read = r
                best_i = i
        if best_olen > 0:
            assembled += best_read[best_olen:]
            reads.pop(best_i)
        else:
           
            assembled += "N" + reads.pop(0)
    return assembled

reconstructed = assemble(samples)
print(f"Reconstructed length: {len(reconstructed)} bases")

out_path = Path("data/reconstructed_sequence.fasta")
with out_path.open("w") as f:
    f.write(">reconstructed\n")
    f.write(reconstructed)
print(f"Saved to {out_path}")

shared = sum(1 for a, b in zip(original_seq, reconstructed) if a == b)
similarity = shared / len(original_seq) * 100
print(f"Similarity to original (approx): {similarity:.2f}%")
