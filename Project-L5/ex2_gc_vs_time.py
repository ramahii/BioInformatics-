from pathlib import Path
import random, time, csv
from Bio import SeqIO

# ---------- helpers (same logic as step 2c) ----------
def gc_content(seq: str) -> float:
    s = seq.upper()
    return 100.0 * (s.count("G") + s.count("C")) / len(s)

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

def generate_reads(sequence: str, n=2000, min_len=100, max_len=150):
    L = len(sequence)
    reads = []
    for _ in range(n):
        k = random.randint(min_len, max_len)
        start = random.randint(0, L - k)
        reads.append(sequence[start:start+k])
    return reads

# ---------- main ----------
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)
seq_path = DATA_DIR / "sequence.fasta"
record = SeqIO.read(seq_path.as_posix(), "fasta")
seq = str(record.seq)

gc = gc_content(seq)
print(f"GC% of current sequence: {gc:.2f}%")

# run multiple times to average timing variability
runs = 3
times = []
for i in range(runs):
    reads = generate_reads(seq)
    t0 = time.perf_counter()
    _ = assemble(reads)
    dt = time.perf_counter() - t0
    times.append(dt)
    print(f"Run {i+1}: {dt:.2f} s")

avg_time = sum(times)/len(times)
print(f"\nAverage reconstruction time: {avg_time:.2f} s")

# log to CSV for later plotting/correlation
out_csv = DATA_DIR / "gc_vs_time.csv"
write_header = not out_csv.exists()
with out_csv.open("a", newline="") as f:
    w = csv.writer(f)
    if write_header:
        w.writerow(["accession", "length", "gc_percent", "avg_time_seconds"])
    w.writerow([record.id, len(seq), f"{gc:.2f}", f"{avg_time:.2f}"])

print(f"Saved to {out_csv}")
