def calculate_nucleotide_percentages(seq):
    seq = seq.upper()
    total_length = len(seq)

    nucleotides = sorted(set(seq))
    
    print("Nucleotide Percentages:")
    for nucleotide in nucleotides:
        count = seq.count(nucleotide)
        percentage = (count / total_length) * 100
        print(f"{nucleotide}: {percentage:.2f}%")

S = "ACGGGCATATGCGC"

calculate_nucleotide_percentages(S)
