"""
Bioinformatics Lab - Exercise 3
File name: ex3.py

This script plots the frequencies of DNA repetitions found
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import sys

def load_dna_sequence(filename="dna_sequence.txt"):
    """Load DNA sequence from file"""
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
    """Find all repetitions in a DNA sequence"""
    print(f"\nSearching for repetitions ({min_length}-{max_length} bp)...")
    repetitions = {}
    
    for length in range(min_length, max_length + 1):
        print(f"  Analyzing {length}-mers...", end='')
        count = 0
        for i in range(len(sequence) - length + 1):
            pattern = sequence[i:i + length]
            if pattern in repetitions:
                repetitions[pattern] += 1
            else:
                repetitions[pattern] = 1
        print(f" Done")
    
    # Keep only patterns that appear more than once
    original_count = len(repetitions)
    repetitions = {pattern: count for pattern, count in repetitions.items() if count > 1}
    print(f"\n✓ Found {len(repetitions)} patterns that repeat (out of {original_count} total)")
    
    return repetitions

def plot_repetition_frequencies(repetitions):
    """
    Create multiple plots to visualize repetition frequencies
    """
    if not repetitions:
        print("✗ No repetitions to plot!")
        return False
    
    print("\nGenerating plots...")
    
    # Sort repetitions by frequency
    sorted_reps = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle('DNA Repetition Frequency Analysis', fontsize=16, fontweight='bold')
    
    # ===== PLOT 1: Top 20 Most Frequent Repetitions =====
    print("  Creating plot 1/4: Top repetitions...")
    ax1 = plt.subplot(2, 2, 1)
    top_20 = sorted_reps[:20]
    patterns = [p[0] for p in top_20]
    frequencies = [p[1] for p in top_20]
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(patterns)))
    bars = ax1.barh(range(len(patterns)), frequencies, color=colors)
    ax1.set_yticks(range(len(patterns)))
    ax1.set_yticklabels(patterns, fontsize=8)
    ax1.set_xlabel('Frequency', fontsize=10)
    ax1.set_title('Top 20 Most Frequent Repetitions', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3)
    
    # ===== PLOT 2: Frequency Distribution by Pattern Length =====
    print("  Creating plot 2/4: Frequencies by length...")
    ax2 = plt.subplot(2, 2, 2)
    by_length = {}
    for pattern, count in repetitions.items():
        length = len(pattern)
        if length not in by_length:
            by_length[length] = []
        by_length[length].append(count)
    
    lengths = sorted(by_length.keys())
    avg_frequencies = [np.mean(by_length[l]) for l in lengths]
    total_frequencies = [sum(by_length[l]) for l in lengths]
    
    x = np.arange(len(lengths))
    width = 0.35
    
    bars1 = ax2.bar(x - width/2, avg_frequencies, width, label='Average Frequency', color='skyblue')
    bars2 = ax2.bar(x + width/2, total_frequencies, width, label='Total Frequency', color='coral')
    
    ax2.set_xlabel('Pattern Length (bp)', fontsize=10)
    ax2.set_ylabel('Frequency', fontsize=10)
    ax2.set_title('Repetition Frequencies by Length', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(lengths)
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # ===== PLOT 3: Distribution of Frequencies (Histogram) =====
    print("  Creating plot 3/4: Frequency distribution...")
    ax3 = plt.subplot(2, 2, 3)
    all_frequencies = list(repetitions.values())
    
    ax3.hist(all_frequencies, bins=30, color='mediumseagreen', edgecolor='black', alpha=0.7)
    ax3.set_xlabel('Frequency', fontsize=10)
    ax3.set_ylabel('Number of Patterns', fontsize=10)
    ax3.set_title('Distribution of Repetition Frequencies', fontsize=12, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # Add statistics text
    stats_text = f'Total patterns: {len(repetitions)}\n'
    stats_text += f'Mean frequency: {np.mean(all_frequencies):.2f}\n'
    stats_text += f'Median frequency: {np.median(all_frequencies):.2f}\n'
    stats_text += f'Max frequency: {max(all_frequencies)}'
    ax3.text(0.65, 0.95, stats_text, transform=ax3.transAxes, 
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ===== PLOT 4: Number of Unique Patterns by Length =====
    print("  Creating plot 4/4: Unique patterns by length...")
    ax4 = plt.subplot(2, 2, 4)
    pattern_counts = [len(by_length[l]) for l in lengths]
    
    bars = ax4.bar(lengths, pattern_counts, color='indianred', edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Pattern Length (bp)', fontsize=10)
    ax4.set_ylabel('Number of Unique Patterns', fontsize=10)
    ax4.set_title('Number of Unique Repetitions by Length', fontsize=12, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    print("\n  Saving plot to file...")
    plt.savefig('repetition_frequencies.png', dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved as 'repetition_frequencies.png'")
    
    plt.close()
    return True

def print_summary_statistics(repetitions):
    """Print summary statistics about repetitions"""
    print(f"\n{'='*60}")
    print(f"SUMMARY STATISTICS")
    print(f"{'='*60}")
    
    frequencies = list(repetitions.values())
    
    print(f"Total unique patterns found: {len(repetitions)}")
    print(f"Total repetition occurrences: {sum(frequencies)}")
    print(f"Average frequency: {np.mean(frequencies):.2f}")
    print(f"Median frequency: {np.median(frequencies):.2f}")
    print(f"Standard deviation: {np.std(frequencies):.2f}")
    print(f"Min frequency: {min(frequencies)}")
    print(f"Max frequency: {max(frequencies)}")
    
    # By length statistics
    by_length = {}
    for pattern, count in repetitions.items():
        length = len(pattern)
        if length not in by_length:
            by_length[length] = []
        by_length[length].append(count)
    
    print(f"\nBy length:")
    for length in sorted(by_length.keys()):
        counts = by_length[length]
        print(f"  {length}-mers: {len(counts)} unique, avg freq: {np.mean(counts):.2f}")

# Main execution
if __name__ == "__main__":
    print("=" * 60)
    print("Exercise 3: Plotting Repetition Frequencies")
    print("=" * 60)
    print()
    
    # Load the DNA sequence
    dna_sequence = load_dna_sequence("dna_sequence.txt")
    
    if dna_sequence:
        # Find repetitions
        repetitions = find_repetitions(dna_sequence, min_length=6, max_length=10)
        
        if repetitions:
            # Print summary statistics
            print_summary_statistics(repetitions)
            
            # Create plots
            success = plot_repetition_frequencies(repetitions)
            
            if success:
                print(f"\n{'='*60}")
                print("Exercise 3 Complete! ✓")
                print("Check 'repetition_frequencies.png' in your folder")
                print("=" * 60)
            else:
                print("\n✗ Failed to create plots")
        else:
            print("✗ No repetitions found to plot")
    else:
        print("\n✗ Cannot proceed without DNA sequence")
        print("Run ex1.py first to generate dna_sequence.txt")
        sys.exit(1)