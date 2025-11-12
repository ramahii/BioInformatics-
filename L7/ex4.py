"""
Bioinformatics Lab - Exercise 4
File name: ex4.py

This script downloads 10 influenza genomes and plots repetition frequencies for each
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez, SeqIO
import time
import os

Entrez.email = "nidal.ramahhii@gmail.com"  

# List of influenza genome accession IDs
INFLUENZA_ACCESSIONS = [
    "CY121680",  # Influenza A virus segment 1
    "CY121681",  # Influenza A virus segment 2
    "CY121682",  # Influenza A virus segment 3
    "CY121683",  # Influenza A virus segment 4
    "CY121684",  # Influenza A virus segment 5
    "CY121685",  # Influenza A virus segment 6
    "CY121686",  # Influenza A virus segment 7
    "CY121687",  # Influenza A virus segment 8
    "CY113677",  # Another influenza strain segment
    "CY113678",  # Another influenza strain segment
]

def fetch_influenza_genome(accession_id, retry=3):
    """
    Fetch an influenza genome from NCBI
    
    Parameters:
    accession_id (str): NCBI accession number
    retry (int): Number of retry attempts
    
    Returns:
    str: DNA sequence or None if failed
    """
    for attempt in range(retry):
        try:
            print(f"  Fetching {accession_id}...", end='')
            handle = Entrez.efetch(db="nucleotide", id=accession_id, 
                                  rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            sequence = str(record.seq)
            print(f" ✓ ({len(sequence)} bp)")
            
            # Be nice to NCBI servers
            time.sleep(0.5)
            
            return sequence
        except Exception as e:
            print(f" ✗ Attempt {attempt + 1} failed")
            if attempt < retry - 1:
                time.sleep(2)
            else:
                print(f"  Error: {e}")
                return None
    return None

def find_repetitions(sequence, min_length=6, max_length=10):
    """Find all repetitions in a DNA sequence"""
    repetitions = {}
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            pattern = sequence[i:i + length]
            if pattern in repetitions:
                repetitions[pattern] += 1
            else:
                repetitions[pattern] = 1
    
    # Keep only patterns that appear more than once
    repetitions = {pattern: count for pattern, count in repetitions.items() if count > 1}
    
    return repetitions

def plot_single_genome_frequencies(repetitions, genome_name, index):
    """
    Create a bar plot for a single genome's repetition frequencies
    
    Returns:
    dict: Statistics about the repetitions
    """
    if not repetitions:
        return None
    
    # Sort and get top 15 repetitions
    sorted_reps = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:15]
    patterns = [p[0] for p in sorted_reps]
    frequencies = [p[1] for p in sorted_reps]
    
    # Calculate statistics
    all_freqs = list(repetitions.values())
    stats = {
        'name': genome_name,
        'total_patterns': len(repetitions),
        'total_occurrences': sum(all_freqs),
        'mean_freq': np.mean(all_freqs),
        'median_freq': np.median(all_freqs),
        'max_freq': max(all_freqs),
        'top_patterns': sorted_reps[:5]
    }
    
    return stats, patterns, frequencies

def plot_all_genomes(genome_data):
    """
    Create comprehensive plots for all genomes
    
    Parameters:
    genome_data: List of tuples (genome_name, repetitions)
    """
    print(f"\nGenerating comprehensive plots for {len(genome_data)} genomes...")
    
    # Filter out failed genomes
    valid_data = [(name, reps) for name, reps in genome_data if reps]
    
    if not valid_data:
        print("✗ No valid genome data to plot!")
        return
    
    n_genomes = len(valid_data)
    
    # ===== FIGURE 1: Individual genome plots (2x5 grid) =====
    print("  Creating individual genome plots...")
    fig1 = plt.figure(figsize=(20, 12))
    fig1.suptitle('Top 15 Repetitions per Influenza Genome', fontsize=16, fontweight='bold')
    
    for idx, (genome_name, repetitions) in enumerate(valid_data, 1):
        ax = plt.subplot(2, 5, idx)
        
        result = plot_single_genome_frequencies(repetitions, genome_name, idx)
        if result:
            stats, patterns, frequencies = result
            
            colors = plt.cm.viridis(np.linspace(0, 1, len(patterns)))
            ax.barh(range(len(patterns)), frequencies, color=colors)
            ax.set_yticks(range(len(patterns)))
            ax.set_yticklabels(patterns, fontsize=7)
            ax.set_xlabel('Frequency', fontsize=8)
            ax.set_title(f'{genome_name}\n({stats["total_patterns"]} patterns)', 
                        fontsize=9, fontweight='bold')
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('influenza_individual_repetitions.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved 'influenza_individual_repetitions.png'")
    plt.close()
    
    # ===== FIGURE 2: Comparative analysis =====
    print("  Creating comparative analysis plots...")
    fig2 = plt.figure(figsize=(16, 10))
    fig2.suptitle('Comparative Repetition Analysis - 10 Influenza Genomes', 
                  fontsize=16, fontweight='bold')
    
    # Collect statistics
    genome_names = []
    total_patterns = []
    mean_freqs = []
    max_freqs = []
    
    for genome_name, repetitions in valid_data:
        if repetitions:
            all_freqs = list(repetitions.values())
            genome_names.append(genome_name)
            total_patterns.append(len(repetitions))
            mean_freqs.append(np.mean(all_freqs))
            max_freqs.append(max(all_freqs))
    
    # Plot 1: Total number of unique patterns per genome
    ax1 = plt.subplot(2, 2, 1)
    colors = plt.cm.Set3(np.linspace(0, 1, len(genome_names)))
    bars = ax1.bar(range(len(genome_names)), total_patterns, color=colors, edgecolor='black')
    ax1.set_xticks(range(len(genome_names)))
    ax1.set_xticklabels(genome_names, rotation=45, ha='right', fontsize=8)
    ax1.set_ylabel('Number of Unique Patterns', fontsize=10)
    ax1.set_title('Total Unique Repetitive Patterns', fontsize=12, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Average frequency per genome
    ax2 = plt.subplot(2, 2, 2)
    bars = ax2.bar(range(len(genome_names)), mean_freqs, color=colors, edgecolor='black')
    ax2.set_xticks(range(len(genome_names)))
    ax2.set_xticklabels(genome_names, rotation=45, ha='right', fontsize=8)
    ax2.set_ylabel('Average Frequency', fontsize=10)
    ax2.set_title('Average Repetition Frequency', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    # Plot 3: Maximum frequency per genome
    ax3 = plt.subplot(2, 2, 3)
    bars = ax3.bar(range(len(genome_names)), max_freqs, color=colors, edgecolor='black')
    ax3.set_xticks(range(len(genome_names)))
    ax3.set_xticklabels(genome_names, rotation=45, ha='right', fontsize=8)
    ax3.set_ylabel('Maximum Frequency', fontsize=10)
    ax3.set_title('Maximum Repetition Frequency', fontsize=12, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # Plot 4: Frequency distribution comparison (box plot)
    ax4 = plt.subplot(2, 2, 4)
    all_freq_data = []
    for genome_name, repetitions in valid_data:
        if repetitions:
            all_freq_data.append(list(repetitions.values()))
    
    bp = ax4.boxplot(all_freq_data, labels=genome_names, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax4.set_xticklabels(genome_names, rotation=45, ha='right', fontsize=8)
    ax4.set_ylabel('Frequency', fontsize=10)
    ax4.set_title('Frequency Distribution Comparison', fontsize=12, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('influenza_comparative_analysis.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved 'influenza_comparative_analysis.png'")
    plt.close()

def save_summary_report(genome_data):
    """Save a text summary report"""
    print("\n  Creating summary report...")
    
    with open("influenza_analysis_report.txt", "w") as f:
        f.write("="*70 + "\n")
        f.write("INFLUENZA GENOME REPETITION ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        for idx, (genome_name, repetitions) in enumerate(genome_data, 1):
            f.write(f"\n{'='*70}\n")
            f.write(f"GENOME {idx}: {genome_name}\n")
            f.write(f"{'='*70}\n")
            
            if repetitions:
                all_freqs = list(repetitions.values())
                f.write(f"Total unique patterns: {len(repetitions)}\n")
                f.write(f"Total repetitions: {sum(all_freqs)}\n")
                f.write(f"Average frequency: {np.mean(all_freqs):.2f}\n")
                f.write(f"Median frequency: {np.median(all_freqs):.2f}\n")
                f.write(f"Max frequency: {max(all_freqs)}\n\n")
                
                f.write("Top 10 most frequent patterns:\n")
                f.write(f"{'Pattern':<15} {'Length':<10} {'Frequency':<12}\n")
                f.write(f"{'-'*15} {'-'*10} {'-'*12}\n")
                
                sorted_reps = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:10]
                for pattern, count in sorted_reps:
                    f.write(f"{pattern:<15} {len(pattern):<10} {count:<12}\n")
            else:
                f.write("Failed to analyze this genome\n")
    
    print("  ✓ Saved 'influenza_analysis_report.txt'")

# Main execution
if __name__ == "__main__":
    print("=" * 70)
    print("Exercise 4: Influenza Genome Repetition Analysis")
    print("=" * 70)
    print()
    
    # Step 1: Fetch genomes
    print(f"Step 1: Downloading {len(INFLUENZA_ACCESSIONS)} influenza genomes from NCBI...")
    print("-" * 70)
    
    genomes = []
    for i, accession in enumerate(INFLUENZA_ACCESSIONS, 1):
        print(f"[{i}/{len(INFLUENZA_ACCESSIONS)}]", end=' ')
        sequence = fetch_influenza_genome(accession)
        
        if sequence:
            genomes.append((accession, sequence))
            # Save each genome
            with open(f"influenza_{accession}.txt", "w") as f:
                f.write(sequence)
        else:
            genomes.append((accession, None))
    
    print(f"\n✓ Successfully downloaded {sum(1 for _, seq in genomes if seq)} genomes")
    
    # Step 2: Analyze repetitions
    print(f"\nStep 2: Analyzing repetitions in each genome...")
    print("-" * 70)
    
    genome_data = []
    for i, (accession, sequence) in enumerate(genomes, 1):
        print(f"[{i}/{len(genomes)}] Analyzing {accession}...", end='')
        
        if sequence:
            repetitions = find_repetitions(sequence, min_length=6, max_length=10)
            genome_data.append((accession, repetitions))
            print(f" ✓ Found {len(repetitions)} unique patterns")
        else:
            genome_data.append((accession, None))
            print(f" ✗ Skipped (no sequence)")
    
    # Step 3: Create plots
    print(f"\nStep 3: Creating visualizations...")
    print("-" * 70)
    plot_all_genomes(genome_data)
    
    # Step 4: Save summary report
    print(f"\nStep 4: Generating summary report...")
    print("-" * 70)
    save_summary_report(genome_data)
    
    print(f"\n{'='*70}")
    print("Exercise 4 Complete! ✓")
    print("="*70)
    print("\nGenerated files:")
    print("  - influenza_individual_repetitions.png")
    print("  - influenza_comparative_analysis.png")
    print("  - influenza_analysis_report.txt")
    print("  - influenza_[ACCESSION].txt (10 genome files)")
    print("="*70)