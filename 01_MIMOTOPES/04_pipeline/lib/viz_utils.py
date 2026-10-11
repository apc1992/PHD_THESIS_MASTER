import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pathlib import Path
from typing import Dict, List

STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

def set_publication_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.titlesize": 11,
        "savefig.dpi": 300,
        "savefig.format": "png",
        "savefig.bbox": "tight"
    })

def generate_restriction_heatmap(
    wt_seq: str, 
    constraints: Dict[int, List[str]], 
    output_path: Path
) -> None:
    set_publication_style()
    
    seq_len = len(wt_seq)
    matrix = np.zeros((len(STANDARD_AAS), seq_len), dtype=bool)

    for pos_idx in range(seq_len):
        allowed = constraints[pos_idx]  # Explicitly defined for all positions now
        for aa in allowed:
            aa_idx = STANDARD_AAS.index(aa)
            matrix[aa_idx, pos_idx] = True

    fig, ax = plt.subplots(figsize=(seq_len * 0.45 + 1.5, 6))
    cmap = ListedColormap(["#FFFFFF", "#F08080"])
    
    ax.imshow(matrix, cmap=cmap, aspect="auto", origin="upper")
    
    ax.set_xticks(np.arange(-0.5, seq_len, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(STANDARD_AAS), 1), minor=True)
    ax.grid(which="minor", color="#D3D3D3", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    xtick_labels = [f"P{i+1}\n({wt_seq[i]})" for i in range(seq_len)]
    ax.set_xticks(np.arange(seq_len))
    ax.set_xticklabels(xtick_labels, fontweight="bold")
    
    ax.set_yticks(np.arange(len(STANDARD_AAS)))
    ax.set_yticklabels(STANDARD_AAS)
    
    ax.set_xlabel("Peptide Position (WT Residue)")
    ax.set_ylabel("Amino Acid Residue")
    ax.set_title(f"PS-CPL Positional Restrictions: {wt_seq}", pad=12, fontweight="bold")

    for pos_idx, wt_aa in enumerate(wt_seq):
        aa_idx = STANDARD_AAS.index(wt_aa)
        ax.plot(pos_idx, aa_idx, marker="o", markersize=4, color="#8B0000")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
