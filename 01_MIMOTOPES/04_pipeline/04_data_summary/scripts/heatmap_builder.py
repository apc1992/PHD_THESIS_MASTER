import logging
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from typing import List

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def generate_donor_heatmaps(
    cytokines_df: pd.DataFrame, 
    wt_seq: str, 
    cyto_cols: List[str], 
    output_dir: Path,
    peptide_name: str
):
    """Generates 4-panel column-wise Z-normalized relative activity heatmaps per donor."""
    output_dir.mkdir(parents=True, exist_ok=True)
    donors = cytokines_df['donor'].unique()
    positions = list(range(1, len(wt_seq) + 1))

    sns.set_theme(style="ticks", context="paper", font_scale=1.2)

    for donor in donors:
        donor_df = cytokines_df[cytokines_df['donor'] == donor]
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True, constrained_layout=True)
        axes = axes.flatten()

        for idx, col in enumerate(cyto_cols[:4]):
            ax = axes[idx]
            matrix = np.full((len(AMINO_ACIDS), len(positions)), np.nan)

            for p_idx, pos in enumerate(positions):
                wt_aa = wt_seq[pos - 1]
                wt_val = donor_df[(donor_df['position'] == pos) & (donor_df['aminoacid'] == wt_aa)][col].values
                base = wt_val[0] if len(wt_val) > 0 and wt_val[0] > 0 else 1e-5

                for a_idx, aa in enumerate(AMINO_ACIDS):
                    val = donor_df[(donor_df['position'] == pos) & (donor_df['aminoacid'] == aa)][col].values
                    if len(val) > 0 and pd.notna(val[0]):
                        matrix[a_idx, p_idx] = val[0] / base

            # Column-wise Z-Score Normalization without zero-variance/empty slice warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                col_means = np.nanmean(matrix, axis=0)
                col_stds = np.nanstd(matrix, axis=0)

            col_means = np.nan_to_num(col_means, nan=0.0)
            col_stds = np.nan_to_num(col_stds, nan=1.0)
            col_stds[col_stds == 0] = 1.0

            z_matrix = np.where(np.isnan(matrix), np.nan, (matrix - col_means) / col_stds)

            heatmap_df = pd.DataFrame(z_matrix, index=AMINO_ACIDS, columns=positions)
            sns.heatmap(
                heatmap_df, cmap="vlag", center=0, ax=ax, 
                cbar_kws={'label': 'Column Z-Score'}, linewidths=0.3, linecolor='gray'
            )

            for pos_idx, wt_aa in enumerate(wt_seq):
                if wt_aa in AMINO_ACIDS:
                    aa_idx = AMINO_ACIDS.index(wt_aa)
                    ax.add_patch(Rectangle((pos_idx, aa_idx), 1, 1, fill=False, edgecolor='black', lw=1.5))

            ax.set_title(f"{col}", fontweight='bold')
            ax.set_xlabel("Peptide Position" if idx >= 2 else "")
            ax.set_ylabel("Amino Acid" if idx % 2 == 0 else "")

        fig.suptitle(f"[{peptide_name}] Donor {donor} Relative Cytokine Activity (Column Z-Scores)", fontsize=16, fontweight='bold')
        out_file = output_dir / f"{peptide_name}_activity_heatmap_donor_{donor}.png"
        plt.savefig(out_file, dpi=300)
        plt.close()
        logging.info(f"Exported heatmap: {out_file}")
