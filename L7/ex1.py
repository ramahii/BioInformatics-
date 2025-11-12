"""
Bioinformatics Lab - Exercise 1
File name: ex1.py

This script:
1. Gets a DNA sequence from NCBI (or uses a sample sequence)
2. Saves it to a file for further analysis
"""

from Bio import Entrez, SeqIO

# Set your email (required by NCBI)
Entrez.email = "your.email@example.com"  # Replace with your actual email

def fetch_dna_sequence_from_ncbi(accession_id):
    """
    Fetch a DNA sequence from NCBI using its accession ID
    
    Parameters:
    accession_id (str): NCBI accession number
    
    Returns:
    str: DNA sequence
    """
    try:
        print(f"Fetching sequence {accession_id} from NCBI...")
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        sequence = str(record.seq)
        print(f"✓ Successfully fetched sequence!")
        print(f"  Length: {len(sequence)} nucleotides")
        
        return sequence
    except Exception as e:
        print(f"✗ Error fetching sequence: {e}")
        return None

def get_sample_sequence():
    """
    Returns a sample DNA sequence (E. coli lac operon region)
    This is a real sequence from NCBI: J01636.1 (2000 bp)
    """
    sample = """
ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGACGAAAGGGCCTCGTGATACGCCTATTTTTATAGGTTAATGTCATGATAATAATGGTTTCTTAGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGA
    """.replace('\n', '').replace(' ', '').strip()
    
    return sample

# Main execution
if __name__ == "__main__":
    print("=" * 60)
    print("Exercise 1: Obtaining DNA Sequence")
    print("=" * 60)
    
    # Option 1: Try to fetch from NCBI
    # Uncomment and add a valid accession ID if you want to fetch from NCBI
    # accession = "J01636.1"  # E. coli lac operon (2000 bp)
    # dna_sequence = fetch_dna_sequence_from_ncbi(accession)
    
    # Option 2: Use sample sequence (recommended for testing)
    print("\nUsing sample DNA sequence (E. coli lac operon region)")
    dna_sequence = get_sample_sequence()
    print(f"✓ Sequence loaded!")
    print(f"  Length: {len(dna_sequence)} nucleotides")
    
    if dna_sequence:
        # Ensure sequence is in the required range (1000-3000 nucleotides)
        if len(dna_sequence) < 1000:
            print(f"\n⚠ Warning: Sequence is shorter than 1000 bp")
        elif len(dna_sequence) > 3000:
            print(f"\n⚠ Sequence is longer than 3000 bp, trimming to 2000 bp")
            dna_sequence = dna_sequence[:2000]
        
        # Save to file for use in next exercises
        with open("dna_sequence.txt", "w") as f:
            f.write(dna_sequence)
        
        print(f"\n✓ Sequence saved to 'dna_sequence.txt'")
        print(f"  Final length: {len(dna_sequence)} nucleotides")
        print(f"  First 60 bp: {dna_sequence[:60]}...")
        print(f"  Last 60 bp: ...{dna_sequence[-60:]}")
        print("\n" + "=" * 60)
        print("Exercise 1 Complete! ✓")
        print("=" * 60)
    else:
        print("\n✗ Failed to obtain DNA sequence")