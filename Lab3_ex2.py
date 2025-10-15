import math

def read_fasta(file_path: str) -> str:
   
    with open(file_path, 'r') as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()


def tm_accurate(dna_seq: str, na_conc: float = 0.05) -> float:
    
    length = len(dna_seq)
    g = dna_seq.count('G')
    c = dna_seq.count('C')
    gc_percent = ((g + c) / length) * 100
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm


def sliding_window_tm(seq: str, window_size: int = 8, na_conc: float = 0.05):
    
    results = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        tm = tm_accurate(window, na_conc)
        results.append((i + 1, window, tm))
    return results


def main():
    
    fasta_path = input("Enter FASTA file path: ").strip()
    na_conc = float(input("Enter Na+ concentration (in M, e.g., 0.05): ") or 0.05)
    window_size = 8  

    
    sequence = read_fasta(fasta_path)
    results = sliding_window_tm(sequence, window_size, na_conc)

    
    print(f"\nDNA sequence (length {len(sequence)}):")
    print(sequence)
    print(f"\nWindow size: {window_size}")
    print(f"Na+ concentration: {na_conc} M\n")
    print("Position\tWindow\t\tTm (°C)")
    print("-" * 40)
    for pos, window, tm in results:
        print(f"{pos:>8}\t{window}\t{tm:6.2f}")


if __name__ == "__main__":
    main()
