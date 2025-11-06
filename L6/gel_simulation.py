import random
import matplotlib.pyplot as plt
import numpy as np

# Generate an arbitrary DNA sequence (length between 1000 and 3000)
sequence_length = random.randint(1000, 3000)
bases = ['A', 'T', 'G', 'C']
dna_sequence = ''.join(random.choice(bases) for _ in range(sequence_length))

# Generate 10 random DNA fragments (between 100–3000 bases)
fragments = []
for _ in range(10):
    frag_length = random.randint(100, min(3000, sequence_length))
    start = random.randint(0, sequence_length - frag_length)
    fragment = dna_sequence[start:start + frag_length]
    fragments.append(fragment)

# Get fragment sizes
fragment_sizes = [len(frag) for frag in fragments]

# Simulate migration: smaller fragments move farther (nonlinear)
min_bp, max_bp = min(fragment_sizes), max(fragment_sizes)
positions = []
for size in fragment_sizes:
    pos = np.interp(np.log(size), [np.log(min_bp), np.log(max_bp)], [0.9, 0.1])
    positions.append(pos)

# Sort fragments by size (larger = top)
sorted_pairs = sorted(zip(fragment_sizes, positions), key=lambda x: x[0], reverse=True)
sorted_sizes, sorted_positions = zip(*sorted_pairs)

# Create the figure
fig, ax = plt.subplots(figsize=(3, 6))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Draw solid black rectangle for the gel background
ax.add_patch(plt.Rectangle((0.1, 0), 0.8, 1, color='black'))

# Draw white border
ax.plot([0.1, 0.9, 0.9, 0.1, 0.1], [0, 0, 1, 1, 0], color='white', linewidth=2)

# Draw bands
for size, pos in zip(sorted_sizes, sorted_positions):
    ax.hlines(pos, 0.2, 0.8, color='white', linewidth=4)
    ax.text(0.82, pos, f"{size} bp", color='white', fontsize=8, va='center')

# Reference labels
ax.text(0.02, 0.9, "3000 bp -", color='black', fontsize=9)
ax.text(0.02, 0.5, "1500 bp -", color='black', fontsize=9)
ax.text(0.02, 0.2, "500 bp -", color='black', fontsize=9)

plt.title("Simulated DNA Gel Electrophoresis", color='black', pad=20)
plt.tight_layout()

# Save and show the figure
plt.savefig("gel_result.png", dpi=300, bbox_inches='tight')
plt.show()
