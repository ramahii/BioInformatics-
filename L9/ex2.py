TEs = {
    "TE1": "ATGCGTACGA",
    "TE2": "TTCGAACTGAAC",
    "TE3": "GGATCCGAT",
    "TE1_copy": "ATGCGTACGA"
}

dna_sequence = """
ACGGGGGGGTGCTTGCATATCGTCCTCTTACATGTCGCTATTCCATGGGAACACGAAATGCGTACGACAGATGCTGGCAATATTTGGTAATAGGACTTCGCGGTCTATGTTCACCCCAGTAAATATAGTAACATGAAGTTCGAACTGAACAGGAAAATAGGCGTAGGATCAGTCCGAAATCTCCTTAGGCAAGCCTAATGGCCATTGAAGCGAAATACTTTCGCTTCCGCCATGCATCGCTGGCGCACCGAAAGTGAATTGCCAAGTTCTTTTCGTATGCGTACGAGCAGGCACCCCCCCGGATCCGATGAGCGGTTCTAAGACGAATTAGTCTATTGTAG
""".replace("\n", "").strip()


def find_transposons(sequence, patterns):
    results = []
    for name, te in patterns.items():
        start = 0
        while True:
            index = sequence.find(te, start)
            if index == -1:
                break
            results.append((name, index, index + len(te)))
            start = index + 1
    return results


if __name__ == "__main__":
    matches = find_transposons(dna_sequence, TEs)
    for name, start, end in matches:
        print(f"{name}: start={start}, end={end}")
    print(f"\nTotal elements detected: {len(matches)}")
