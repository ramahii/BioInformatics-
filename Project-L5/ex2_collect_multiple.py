from pathlib import Path
import random, time, csv
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"  # ← use your email!

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)
CSV_PATH = DATA_DIR / "gc_vs_time.csv"

# ---------- helpers ----------
def fetch_random_sequence(min_len=1000, max_len=3000):
    query = f"(biomol_genomic[PROP] OR biomol_mrna[PROP]) AND {min_len}:{max_len}[SLEN]"
    with Entrez.esearch(db="nucleotide", term=query, retmax=200) as h:
        ids = Entrez.read(h)["IdList"]
    gi = random.choice(ids)
    with Entrez.efetch(db="nucleotide", id=gi, rettype="fasta", retmode="text") as h:
        return SeqIO.read(h, "fasta")

def gc_content(seq):
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)

def overlap(a, b, min_overlap=20):
    start = 0
    while True:
        start = a.find(b[:min_overlap], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def assemble(reads):
    reads = reads.copy()
    assembled = reads.pop(0)
    while reads:
        best_olen, best_i = 0, None
        for i, r in enumerate(reads):
            olen = overlap(assembled, r)
            if olen > best_olen:
                best_olen, best_i = olen, i
        if best_olen > 0:
            assembled += reads[best_i][best_olen:]
            reads.pop(best_i)
        else:
            assembled += "N" + reads.pop(0)
    return assembled

def generate_reads(seq, n=2000, min_len=100, max_len=150):
    L = len(seq)
    reads = []
    for _ in range(n):
        k = random.randint(min_len, max_len)
        start = random.randint(0, L - k)
        reads.append(seq[start:start + k])
    return reads

# ---------- experiment ----------
def measure_one():
    rec = fetch_random_sequence()
    seq = str(rec.seq)
    gc = gc_content(seq)

    reads = generate_reads(seq)
    t0 = time.perf_counter()
    assemble(reads)
    dt = time.perf_counter() - t0
    return rec.id, len(seq), gc, dt

# ---------- run multiple ----------
N = 5  # number of random sequences to test
with CSV_PATH.open("a", newline="") as f:
    w = csv.writer(f)
    if f.tell() == 0:
        w.writerow(["accession", "length", "gc_percent", "avg_time_seconds"])

    for i in range(N):
        acc, L, gc, t = measure_one()
        w.writerow([acc, L, f"{gc:.2f}", f"{t:.2f}"])
        print(f"{i+1}/{N}: {acc} | len={L} | GC={gc:.2f}% | time={t:.2f}s")

print(f"\n✅ Added {N} new sequences to {CSV_PATH}")
