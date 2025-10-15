import math

def tm_basic(seq):
   
    seq = seq.upper()
    A = seq.count('A')
    T = seq.count('T')
    G = seq.count('G')
    C = seq.count('C')
    
    tm = 4 * (G + C) + 2 * (A + T)
    return tm

def tm_advanced(seq, Na_conc=0.05):
    
    seq = seq.upper()
    length = len(seq)
    G = seq.count('G')
    C = seq.count('C')
    GC_percent = ((G + C) / length) * 100

    tm = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * GC_percent - (600 / length)
    return tm

if __name__ == "__main__":
    dna_seq = input("Enter a DNA sequence: ").strip()
    if not dna_seq:
        print("Please enter a valid DNA sequence!")
        exit()

    basic_tm = tm_basic(dna_seq)
    adv_tm = tm_advanced(dna_seq)

    print("\n--- DNA Melting Temperature Results ---")
    print(f"DNA Sequence: {dna_seq}")
    print(f"Length: {len(dna_seq)} bases")
    print(f"Basic Formula (Wallace rule): {basic_tm:.2f} °C")
    print(f"Advanced Formula: {adv_tm:.2f} °C")
