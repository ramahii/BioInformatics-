import random

# --- Step 1: Define transposable elements ---
TEs = {
    "TE1": "ATGCGTACGA",
    "TE2": "TTCGAACTGAAC",
    "TE3": "GGATCCGAT"
}

def random_dna(length):
    """Generates a random DNA sequence of given length."""
    return "".join(random.choice("ATGC") for _ in range(length))

def create_sequence(base_length=300):
    """Creates a random DNA sequence and inserts transposable elements."""
    base_sequence = random_dna(base_length)
    sequence_list = list(base_sequence)
    positions = []

    # Insert each TE once
    for name, te_seq in TEs.items():
        pos = random.randint(0, len(sequence_list))
        sequence_list[pos:pos] = te_seq
        positions.append((name, pos, pos + len(te_seq)))

    # Optional: insert TE1 again for a 4th TE
    name = "TE1_copy"
    te_seq = TEs["TE1"]
    pos = random.randint(0, len(sequence_list))
    sequence_list[pos:pos] = te_seq
    positions.append((name, pos, pos + len(te_seq)))

    final_sequence = "".join(sequence_list)
    return final_sequence, positions


# ----- RUNNING THE EXERCISE -----
dna, positions = create_sequence()

print("Final sequence length:", len(dna))
print("\nTransposable element positions:")

for name, start, end in positions:
    print(f"{name}: start={start}, end={end}")

print("\nDNA Sequence:")
print(dna)
