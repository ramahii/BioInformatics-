from collections import Counter

def read_fasta(filename):
    """Reads a FASTA file and returns the concatenated sequence as a string."""
    sequence = []
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):  
                sequence.append(line.strip().upper())  
    return "".join(sequence)

def find_alphabet(sequence):
    """Finds the unique symbols (alphabet) in the sequence."""
    return sorted(set(sequence))

def find_percentages(sequence):
    """Finds relative percentages of each symbol in the sequence."""
    length = len(sequence)
    counts = Counter(sequence)
    percentages = {symbol: (count / length) * 100 for symbol, count in counts.items()}
    return percentages

if __name__ == "__main__":
    fasta_file = "example.fasta"  
    seq = read_fasta(fasta_file)

    alphabet = find_alphabet(seq)
    percentages = find_percentages(seq)

    print("Alphabet:", alphabet)
    print("Percentages:")
    for symbol, perc in percentages.items():
        print(f"  {symbol}: {perc:.2f}%")