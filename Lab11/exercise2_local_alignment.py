import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
import numpy as np
import os
import urllib.request
import threading

# Default sequences for testing (small example)
DEFAULT_S1 = "GATTACAGATTACA"
DEFAULT_S2 = "GCATGCTAGATTACA"

DEFAULT_MATCH = 2
DEFAULT_MISMATCH = -1
DEFAULT_GAP = -1

CHUNK_SIZE = 500  # Process genomes in 500 bp chunks


def load_fasta_file(filepath):
    """Load FASTA file and return sequence as string"""
    sequence = ""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('>'):
                    sequence += line.upper()
        return sequence
    except Exception as e:
        raise Exception(f"Error loading file: {str(e)}")


def download_fasta_from_ncbi(accession_id, filepath, callback=None):
    """
    Download FASTA sequence from NCBI using Entrez Direct API
    accession_id: NCBI accession number (e.g., NC_045512.2 for COVID-19)
    filepath: where to save the file
    callback: function to call with progress updates
    """
    try:
        url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={accession_id}"
        
        if callback:
            callback(f"Downloading {accession_id}...")
        
        urllib.request.urlretrieve(url, filepath)
        
        if callback:
            callback(f"Successfully downloaded to {filepath}")
        
        return True
    except Exception as e:
        if callback:
            callback(f"Error downloading: {str(e)}")
        raise Exception(f"Failed to download from NCBI: {str(e)}")


def smith_waterman_chunk(seq1, seq2, match_score, mismatch_score, gap_penalty):
    """Standard Smith-Waterman for local alignment (small sequences)"""
    n = len(seq1)
    m = len(seq2)

    # Score and traceback matrices
    score = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[None] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    # Fill matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag_score = score[i - 1][j - 1] + match_score
            else:
                diag_score = score[i - 1][j - 1] + mismatch_score

            up_score = score[i - 1][j] + gap_penalty
            left_score = score[i][j - 1] + gap_penalty

            max_val = max(0, diag_score, up_score, left_score)
            score[i][j] = max_val

            if max_val > max_score:
                max_score = max_val
                max_i, max_j = i, j

            if max_val == diag_score and max_val > 0:
                traceback[i][j] = "diag"
            elif max_val == up_score and max_val > 0:
                traceback[i][j] = "up"
            elif max_val == left_score and max_val > 0:
                traceback[i][j] = "left"

    # Traceback from max_score position
    aligned1 = []
    aligned2 = []
    i, j = max_i, max_j
    path = [(i, j)]

    while i > 0 or j > 0:
        direction = traceback[i][j] if i > 0 and j > 0 else None
        if direction == "diag":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        elif direction == "left":
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            break
        path.append((i, j))

    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    path = list(reversed(path))

    return aligned1, aligned2, score, path, max_score


def chunked_local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, chunk_size=CHUNK_SIZE):
    """
    Divide large sequences into overlapping chunks and perform local alignment.
    Returns aggregated results showing region similarities.
    """
    all_alignments = []
    all_scores = []
    
    # Process seq1 in chunks
    for i in range(0, len(seq1), chunk_size):
        chunk1 = seq1[i:i + chunk_size]
        chunk_start1 = i
        chunk_end1 = i + len(chunk1)
        
        # Process seq2 in chunks
        for j in range(0, len(seq2), chunk_size):
            chunk2 = seq2[j:j + chunk_size]
            chunk_start2 = j
            chunk_end2 = j + len(chunk2)
            
            # Perform local alignment on chunk pair
            try:
                aligned1, aligned2, score_matrix, path, max_score = smith_waterman_chunk(
                    chunk1, chunk2, match_score, mismatch_score, gap_penalty
                )
                
                if max_score > 0:  # Only record if alignment found
                    all_alignments.append({
                        'seq1': aligned1,
                        'seq2': aligned2,
                        'score': max_score,
                        'start1': chunk_start1,
                        'end1': chunk_end1,
                        'start2': chunk_start2,
                        'end2': chunk_end2,
                        'matrix': score_matrix,
                        'path': path
                    })
                    all_scores.append(max_score)
            except:
                continue
    
    return all_alignments, all_scores


class SmithWatermanApp:
    def __init__(self, master):
        self.master = master
        master.title("Smith-Waterman Local Alignment - Large Genomes")
        master.geometry("1400x900")

        # ====== Top frame: File loading and parameters ======
        top_frame = ttk.Frame(master, padding=5)
        top_frame.pack(side=tk.TOP, fill=tk.X)

        # File loading section
        file_frame = ttk.LabelFrame(top_frame, text="Load FASTA Files or Use Defaults")
        file_frame.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        ttk.Button(file_frame, text="Load Seq 1 (FASTA)", command=self.load_seq1).pack(side=tk.LEFT, padx=5)
        ttk.Button(file_frame, text="Load Seq 2 (FASTA)", command=self.load_seq2).pack(side=tk.LEFT, padx=5)
        ttk.Label(file_frame, text="Or use default sequences").pack(side=tk.LEFT, padx=10)

        # Sequences display
        seq_frame = ttk.LabelFrame(top_frame, text="Sequences Info")
        seq_frame.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        self.seq1_label = ttk.Label(seq_frame, text=f"Seq 1: {len(DEFAULT_S1)} bp")
        self.seq1_label.pack(anchor="w")
        self.seq2_label = ttk.Label(seq_frame, text=f"Seq 2: {len(DEFAULT_S2)} bp")
        self.seq2_label.pack(anchor="w")

        # Parameters
        params_frame = ttk.LabelFrame(top_frame, text="Parameters")
        params_frame.grid(row=0, column=2, padx=15, pady=5, sticky="n")

        ttk.Label(params_frame, text="Match:").grid(row=0, column=0, sticky="e")
        self.e_match = ttk.Entry(params_frame, width=8)
        self.e_match.grid(row=0, column=1, padx=3, pady=2)
        self.e_match.insert(0, str(DEFAULT_MATCH))

        ttk.Label(params_frame, text="Mismatch:").grid(row=1, column=0, sticky="e")
        self.e_mismatch = ttk.Entry(params_frame, width=8)
        self.e_mismatch.grid(row=1, column=1, padx=3, pady=2)
        self.e_mismatch.insert(0, str(DEFAULT_MISMATCH))

        ttk.Label(params_frame, text="Gap:").grid(row=2, column=0, sticky="e")
        self.e_gap = ttk.Entry(params_frame, width=8)
        self.e_gap.grid(row=2, column=1, padx=3, pady=2)
        self.e_gap.insert(0, str(DEFAULT_GAP))

        ttk.Label(params_frame, text="Chunk Size:").grid(row=3, column=0, sticky="e")
        self.e_chunk = ttk.Entry(params_frame, width=8)
        self.e_chunk.grid(row=3, column=1, padx=3, pady=2)
        self.e_chunk.insert(0, str(CHUNK_SIZE))

        # Align button
        self.align_button = ttk.Button(top_frame, text="Align Genomes", command=self.run_alignment)
        self.align_button.grid(row=0, column=3, padx=10, pady=5, sticky="n")

        # ====== Middle frame: plots and text ======
        middle_frame = ttk.Frame(master, padding=5)
        middle_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Left: Heatmap visualization
        self.fig_heatmap = Figure(figsize=(6, 5), dpi=100)
        self.ax_heatmap = self.fig_heatmap.add_subplot(111)
        self.canvas_heatmap = FigureCanvasTkAgg(self.fig_heatmap, master=middle_frame)
        self.canvas_heatmap_widget = self.canvas_heatmap.get_tk_widget()
        self.canvas_heatmap_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        # Right: Text results
        right_frame = ttk.Frame(middle_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5)

        ttk.Label(right_frame, text="Alignment Results:").pack(anchor="w")
        self.text = tk.Text(right_frame, width=60, height=30, font=('Courier', 8))
        self.text.pack(fill=tk.BOTH, expand=True)

        # Store sequences
        self.seq1 = DEFAULT_S1
        self.seq2 = DEFAULT_S2
        self.alignments = None

        # Initial run
        self.run_alignment()

    def load_seq1(self):
        filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa *.fasta"), ("All files", "*.*")])
        if filepath:
            try:
                self.seq1 = load_fasta_file(filepath)
                self.seq1_label.config(text=f"Seq 1: {len(self.seq1)} bp (from {os.path.basename(filepath)})")
                messagebox.showinfo("Success", f"Loaded {len(self.seq1)} bp")
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def load_seq2(self):
        filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa *.fasta"), ("All files", "*.*")])
        if filepath:
            try:
                self.seq2 = load_fasta_file(filepath)
                self.seq2_label.config(text=f"Seq 2: {len(self.seq2)} bp (from {os.path.basename(filepath)})")
                messagebox.showinfo("Success", f"Loaded {len(self.seq2)} bp")
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def download_covid19(self):
        """Download COVID-19 genome from NCBI"""
        self.status_label.config(text="Downloading COVID-19...", foreground="blue")
        self.master.update()
        
        def download_thread():
            try:
                filepath = "covid19.fasta"
                download_fasta_from_ncbi("NC_045512.2", filepath, self.update_status)
                self.seq1 = load_fasta_file(filepath)
                self.seq1_label.config(text=f"Seq 1: {len(self.seq1)} bp (COVID-19)")
                self.status_label.config(text="COVID-19 loaded successfully!", foreground="green")
                messagebox.showinfo("Success", f"Downloaded and loaded COVID-19: {len(self.seq1)} bp")
            except Exception as e:
                self.status_label.config(text="Error downloading", foreground="red")
                messagebox.showerror("Error", f"Failed to download COVID-19: {str(e)}")
        
        thread = threading.Thread(target=download_thread, daemon=True)
        thread.start()

    def download_influenza(self):
        """Download Influenza genome from NCBI"""
        self.status_label.config(text="Downloading Influenza...", foreground="blue")
        self.master.update()
        
        def download_thread():
            try:
                filepath = "influenza.fasta"
                # Using a common influenza A strain (H1N1)
                download_fasta_from_ncbi("NC_002016.1", filepath, self.update_status)
                self.seq2 = load_fasta_file(filepath)
                self.seq2_label.config(text=f"Seq 2: {len(self.seq2)} bp (Influenza)")
                self.status_label.config(text="Influenza loaded successfully!", foreground="green")
                messagebox.showinfo("Success", f"Downloaded and loaded Influenza: {len(self.seq2)} bp")
            except Exception as e:
                self.status_label.config(text="Error downloading", foreground="red")
                messagebox.showerror("Error", f"Failed to download Influenza: {str(e)}")
        
        thread = threading.Thread(target=download_thread, daemon=True)
        thread.start()

    def update_status(self, message):
        """Update status label from download thread"""
        self.status_label.config(text=message, foreground="blue")
        self.master.update()

    def get_params(self):
        try:
            match = int(self.e_match.get())
        except ValueError:
            match = DEFAULT_MATCH

        try:
            mismatch = int(self.e_mismatch.get())
        except ValueError:
            mismatch = DEFAULT_MISMATCH

        try:
            gap = int(self.e_gap.get())
        except ValueError:
            gap = DEFAULT_GAP

        try:
            chunk = int(self.e_chunk.get())
        except ValueError:
            chunk = CHUNK_SIZE

        return match, mismatch, gap, chunk

    def run_alignment(self):
        match, mismatch, gap, chunk_size = self.get_params()

        # Perform chunked local alignment
        alignments, scores = chunked_local_alignment(
            self.seq1, self.seq2, match, mismatch, gap, chunk_size
        )

        self.alignments = alignments

        # Update text results
        self.text.delete("1.0", tk.END)

        self.text.insert(tk.END, f"=== LOCAL ALIGNMENT RESULTS ===\n")
        self.text.insert(tk.END, f"Seq 1 length: {len(self.seq1)} bp\n")
        self.text.insert(tk.END, f"Seq 2 length: {len(self.seq2)} bp\n")
        self.text.insert(tk.END, f"Chunk size: {chunk_size} bp\n\n")

        self.text.insert(tk.END, f"Total high-scoring regions found: {len(alignments)}\n\n")

        if alignments:
            total_score = sum(scores)
            avg_score = total_score / len(scores) if scores else 0
            self.text.insert(tk.END, f"Total alignment score: {total_score}\n")
            self.text.insert(tk.END, f"Average score per region: {avg_score:.1f}\n")
            self.text.insert(tk.END, f"Max score: {max(scores)}\n\n")

            # Show top 5 alignments
            self.text.insert(tk.END, "=== TOP 5 LOCAL ALIGNMENTS ===\n\n")
            sorted_alignments = sorted(alignments, key=lambda x: x['score'], reverse=True)[:5]

            for idx, aln in enumerate(sorted_alignments, 1):
                self.text.insert(tk.END, f"Region {idx}:\n")
                self.text.insert(tk.END, f"  Score: {aln['score']}\n")
                self.text.insert(tk.END, f"  Seq1 region: {aln['start1']}-{aln['end1']} bp\n")
                self.text.insert(tk.END, f"  Seq2 region: {aln['start2']}-{aln['end2']} bp\n")
                self.text.insert(tk.END, f"  Alignment:\n")
                self.text.insert(tk.END, f"    {aln['seq1'][:50]}{'...' if len(aln['seq1']) > 50 else ''}\n")
                self.text.insert(tk.END, f"    {aln['seq2'][:50]}{'...' if len(aln['seq2']) > 50 else ''}\n\n")
        else:
            self.text.insert(tk.END, "No significant local alignments found.\n")

        # Update heatmap visualization
        self.update_heatmap()

    def update_heatmap(self):
        """Visualize alignment scores across genome regions"""
        self.ax_heatmap.clear()

        if not self.alignments:
            self.ax_heatmap.text(0.5, 0.5, 'No alignments found', 
                               ha='center', va='center', transform=self.ax_heatmap.transAxes)
            self.fig_heatmap.tight_layout()
            self.canvas_heatmap.draw_idle()
            return

        # Create heatmap data: genome regions vs alignment scores
        seq1_len = len(self.seq1)
        seq2_len = len(self.seq2)
        grid_size = 50  # Create 50x50 grid

        heatmap = np.zeros((grid_size, grid_size))

        for aln in self.alignments:
            # Map alignment positions to grid
            i1 = int((aln['start1'] / seq1_len) * grid_size)
            j1 = int((aln['start2'] / seq2_len) * grid_size)
            i2 = min(int((aln['end1'] / seq1_len) * grid_size), grid_size - 1)
            j2 = min(int((aln['end2'] / seq2_len) * grid_size), grid_size - 1)

            # Normalize score to 0-1 range
            norm_score = min(aln['score'] / 100, 1.0)
            heatmap[i1:i2+1, j1:j2+1] += norm_score

        im = self.ax_heatmap.imshow(heatmap, cmap='YlOrRd', origin='upper', aspect='auto')
        self.ax_heatmap.set_xlabel("Sequence 2 Regions")
        self.ax_heatmap.set_ylabel("Sequence 1 Regions")
        self.ax_heatmap.set_title("Genome Similarity Heatmap (Local Alignment Regions)")

        import matplotlib.pyplot as plt
        plt.colorbar(im, ax=self.ax_heatmap, label="Alignment Score Density")

        self.fig_heatmap.tight_layout()
        self.canvas_heatmap.draw_idle()


if __name__ == "__main__":
    root = tk.Tk()
    app = SmithWatermanApp(root)
    root.mainloop()