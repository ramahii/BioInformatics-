import sys
import collections
import matplotlib.pyplot as plt

genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def load_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
    return ''.join(seq)

def codon_counts(seq):
    counts = collections.Counter()
    seq = seq.replace('U','T')
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if len(codon)==3:
            counts[codon] += 1
    return counts

def codon_freq(counts):
    total = sum(counts.values())
    return {codon: counts[codon]/total for codon in counts}

def top_n(freq_dict, n=10):
    return sorted(freq_dict.items(), key=lambda x: x[1], reverse=True)[:n]

def aa_counts_from_codon_counts(counts):
    aa_counts = collections.Counter()
    for codon, cnt in counts.items():
        aa = genetic_code.get(codon, None)
        if aa and aa!='Stop':
            aa_counts[aa] += cnt
    return aa_counts

def plot_top(codons, freqs, title, fname):
    plt.figure(figsize=(10,6))
    plt.bar(codons, freqs)
    plt.title(title)
    plt.ylabel('Frequency')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

def main(covid_fasta, flu_fasta):
    covid_seq = load_fasta(covid_fasta)
    flu_seq = load_fasta(flu_fasta)

    covid_counts = codon_counts(covid_seq)
    flu_counts   = codon_counts(flu_seq)

    covid_freq   = codon_freq(covid_counts)
    flu_freq     = codon_freq(flu_counts)

    covid_top10  = top_n(covid_freq, 10)
    flu_top10    = top_n(flu_freq, 10)

    covid_codons, covid_freqs = zip(*covid_top10)
    flu_codons,   flu_freqs   = zip(*flu_top10)

    plot_top(covid_codons, covid_freqs, 'Top 10 Codons – COVID-19 genome', 'covid_top10.png')
    plot_top(flu_codons,   flu_freqs,   'Top 10 Codons – Influenza genome', 'flu_top10.png')

    print("Top 10 codons COVID-19:", covid_top10)
    print("Top 10 codons Influenza:",   flu_top10)

    # comparison part c:
    common = set(codon for codon, _ in covid_top10) & set(codon for codon, _ in flu_top10)
    print("Codons common in top10 of both:", common)

    # part d – top3 amino acids each
    covid_aa_counts = aa_counts_from_codon_counts(covid_counts)
    flu_aa_counts   = aa_counts_from_codon_counts(flu_counts)

    covid_top3_aa = covid_aa_counts.most_common(3)
    flu_top3_aa   = flu_aa_counts.most_common(3)

    print("Top 3 amino acids COVID-19:", covid_top3_aa)
    print("Top 3 amino acids Influenza:",   flu_top3_aa)

    # part e – prompt formulation
    covid_aa3 = [x for x, _ in covid_top3_aa]
    flu_aa3   = [x for x, _ in flu_top3_aa]
    prompt = (f"Which foods are rich in the amino acids {', '.join(covid_aa3+flu_aa3)}? "
              f"Please list food sources and approximate amounts per 100 g.")
    print("\nPrompt to ask AI:")
    print(prompt)

if __name__=='__main__':
    if len(sys.argv)!=3:
        print("Usage: python script.py covid_genome.fasta influenza_genome.fasta")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
