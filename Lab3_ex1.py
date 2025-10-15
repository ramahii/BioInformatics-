import math

def tm_simple(dna_seq: str) -> float:
   
    dna_seq = dna_seq.upper()
    g = dna_seq.count('G')
    c = dna_seq.count('C')
    a = dna_seq.count('A')
    t = dna_seq.count('T')
    tm = 4 * (g + c) + 2 * (a + t)
    return tm


def tm_accurate(dna_seq: str, na_conc: float = 0.05) -> float:
  
    dna_seq = dna_seq.upper()
    length = len(dna_seq)
    g = dna_seq.count('G')
    c = dna_seq.count('C')
    gc_percent = ((g + c) / length) * 100
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm



dna = input("Enter DNA sequence: ").strip().upper()
choice = input("Choose formula (1 for simple, 2 for accurate): ")

if choice == '1':
    result = tm_simple(dna)
    print(f"Melting Temperature (Tm): {result:.2f} °C (Simple Formula)")
elif choice == '2':
    na_conc = float(input("Enter Na+ concentration (in M, e.g., 0.05): "))
    result = tm_accurate(dna, na_conc)
    print(f"Melting Temperature (Tm): {result:.2f} °C (Accurate Formula)")
else:
    print("Invalid choice.")
