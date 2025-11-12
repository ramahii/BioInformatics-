"""
Bioinformatics Lab - Exercise 2
File name: ex2.py

This script detects repetitions (between 6 and 10 nucleotides) in a DNA sequence
"""

def load_dna_sequence(filename="dna_sequence.txt"):
    """
    Load DNA sequence from file
    
    Parameters:
    filename (str): Name of the file containing the DNA sequence
    
    Returns:
    str: DNA sequence
    """
    try:
        with open(filename, 'r') as f:
            sequence = f.read().strip()
        print(f"✓ Loaded sequence from '{filename}'")
        print(f"  Length: {len(sequence)} nucleotides")
        return sequence
    except FileNotFoundError:
        print(f"✗ Error: File '{filename}' not found!")
        print("  Please run ex1.py first to generate the sequence file.")
        return None

def find_repetitions(sequence, min_length=6, max_length=10):
    """
    Find all repetitions in a DNA sequence between min_length and max_length
    
    Parameters:
    sequence (str): DNA sequence
    min_length (int): Minimum length of repetitions to find
    max_length (int): Maximum length of repetitions to find
    
    Returns:
    dict: Dictionary with patterns as keys and their frequencies as values
    """
    repetitions = {}
    
    # For each possible length of repetition
    for length in range(min_length, max_length + 1):
        print(f"  Searching for {length}-mers...")
        
        # Slide through the sequence
        for i in range(len(sequence) - length + 1):
            pattern = sequence[i:i + length]
            
            # Count this pattern
            if pattern in repetitions:
                repetitions[pattern] += 1
            else:
                repetitions[pattern] = 1
    
    # Filter to keep only patterns that appear more than once (actual repetitions)
    repetitions = {pattern: count for pattern, count in repetitions.items() if count > 1}
    
    return repetitions

def analyze_repetitions(repetitions):
    """
    Analyze and display statistics about found repetitions
    
    Parameters:
    repetitions (dict): Dictionary of patterns and their frequencies
    """
    if not repetitions:
        print("\n✗ No repetitions found!")
        return
    
    print(f"\n{'='*60}")
    print(f"REPETITION ANALYSIS")
    print(f"{'='*60}")
    print(f"Total unique patterns found: {len(repetitions)}")
    print(f"Total repetitions: {sum(repetitions.values())}")
    
    # Group by length
    by_length = {}
    for pattern, count in repetitions.items():
        length = len(pattern)
        if length not in by_length:
            by_length[length] = []
        by_length[length].append((pattern, count))
    
    print(f"\nRepetitions by length:")
    for length in sorted(by_length.keys()):
        patterns = by_length[length]
        total_count = sum(count for _, count in patterns)
        print(f"  {length}-mers: {len(patterns)} unique patterns, {total_count} total occurrences")
    
    # Show top 10 most frequent repetitions
    print(f"\n{'='*60}")
    print(f"TOP 10 MOST FREQUENT REPETITIONS")
    print(f"{'='*60}")
    sorted_reps = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:10]
    
    print(f"{'Pattern':<15} {'Length':<10} {'Frequency':<12}")
    print(f"{'-'*15} {'-'*10} {'-'*12}")
    for pattern, count in sorted_reps:
        print(f"{pattern:<15} {len(pattern):<10} {count:<12}")
    
    return repetitions

# Main execution
if __name__ == "__main__":
    print("=" * 60)
    print("Exercise 2: Detecting DNA Repetitions")
    print("=" * 60)
    print()
    
    # Load the DNA sequence
    dna_sequence = load_dna_sequence("dna_sequence.txt")
    
    if dna_sequence:
        print(f"\nSearching for repetitions (6-10 nucleotides)...")
        
        # Find repetitions
        repetitions = find_repetitions(dna_sequence, min_length=6, max_length=10)
        
        # Analyze and display results
        analyze_repetitions(repetitions)
        
        # Save results to file
        with open("repetitions.txt", "w") as f:
            f.write("DNA Repetition Analysis Results\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total unique patterns: {len(repetitions)}\n")
            f.write(f"Total repetitions: {sum(repetitions.values())}\n\n")
            f.write("Pattern\tLength\tFrequency\n")
            f.write("-" * 40 + "\n")
            sorted_reps = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)
            for pattern, count in sorted_reps:
                f.write(f"{pattern}\t{len(pattern)}\t{count}\n")
        
        print(f"\n{'='*60}")
        print(f"✓ Results saved to 'repetitions.txt'")
        print(f"{'='*60}")
        print("\nExercise 2 Complete! ✓")
        print("=" * 60)
    else:
        print("\n✗ Cannot proceed without DNA sequence")