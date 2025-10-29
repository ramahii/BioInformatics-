from pathlib import Path
import random
from Bio import Entrez, SeqIO


Entrez.email = "nidalwadii@gmail.com"  

OUT_DIR = Path("data")
OUT_DIR.mkdir(exist_ok=True)
OUT_FASTA = OUT_DIR / "sequence.fasta"

def fetch_random_sequence(min_len=1000, max_len=3000, term="(biomol_genomic[PROP] OR biomol_mrna[PROP])"):
   
    query = f"{term} AND {min_len}:{max_len}[SLEN]"
    with Entrez.esearch(db="nucleotide", term=query, retmax=200) as handle:
        res = Entrez.read(handle)
    ids = res.get("IdList", [])
    if not ids:
        raise RuntimeError("No sequences found.")
    gi = random.choice(ids)
    with Entrez.efetch(db="nucleotide", id=gi, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    return record

def main():
    rec = fetch_random_sequence()
    SeqIO.write(rec, OUT_FASTA.as_posix(), "fasta")
    print(f"Saved to {OUT_FASTA}")
    print(f"> {rec.id}")
    print(f"Length: {len(rec.seq)} bases")
    print(str(rec.seq)[:120] + ("..." if len(rec.seq) > 120 else ""))

if __name__ == "__main__":
    main()
