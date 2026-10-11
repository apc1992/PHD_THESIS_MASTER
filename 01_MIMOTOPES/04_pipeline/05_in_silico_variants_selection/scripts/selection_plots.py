import logging
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

logger = logging.getLogger(__name__)

# Manuscript-Consistent Category Colors
CATEGORY_COLORS = {
    '1AA': '#1f77b4',     # Navy Blue
    '2AA': '#e6550d',     # Coral / Orange
    'PS-CPL': '#756bb1'   # Purple
}
BACKGROUND_COLOR = '#e0e0e0'

def setup_publication_style():
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({
        'font.sans-serif': 'Helvetica',
        'font.family': 'sans-serif',
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8.5,
        'figure.titlesize': 13
    })

def plot_selection_dashboard(df: pd.DataFrame, selected_df: pd.DataFrame, out_path: Path):
    """Generates an optimized publication dashboard plotting Affinity vs %EL Rank for Panels A-C 
    and a log-transformed 3D profile for Panel D, with SurVaxM sequence highlighted in Panel A.
    """
    setup_publication_style()
    fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    
    el_col = 'HLA-A*02:01_%Rank_EL' if 'HLA-A*02:01_%Rank_EL' in df.columns else 'EL_rank_A0201'
    aff_col = 'HLA-A*02:01_Aff(nM)' if 'HLA-A*02:01_Aff(nM)' in df.columns else 'BA_nM_A0201'
    
    selected_mimotopes = set(selected_df['mimotope'])
    
    # Helper function for 2D Affinity vs %EL Rank panels
    def plot_arm_panel(ax, category_key, title_text, color):
        arm_df = df[df['category'].str.upper().str.replace('-', '') == category_key.replace('-', '')].copy()
        if arm_df.empty:
            return
        
        arm_df['is_selected'] = arm_df['mimotope'].isin(selected_mimotopes)
        bg_df = arm_df[~arm_df['is_selected']]
        sel_df = arm_df[arm_df['is_selected']]
        
        # Plot background unselected candidates
        ax.scatter(bg_df[el_col], bg_df[aff_col], c=BACKGROUND_COLOR, s=18, alpha=0.35, label='Unselected')
        
        # Plot selected candidates
        ax.scatter(sel_df[el_col], sel_df[aff_col], c=color, s=45, alpha=0.9, edgecolor='k', linewidth=0.5, label='Selected Top 50')
        
        # Annotate SurVaxM mimotope (QMFFCFKEL) in 1-AA Panel
        survax_match = arm_df[arm_df['mimotope'].str.upper() == 'QMFFCFKEL']
        if not survax_match.empty:
            sx = survax_match[el_col].values[0]
            sy = survax_match[aff_col].values[0]
            ax.annotate(
                'QMFFCFKEL',
                xy=(sx, sy),
                xytext=(20, 25),
                textcoords='offset points',
                fontsize=8.5,
                fontweight='bold',
                color='#000000',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffffff', edgecolor=color, linewidth=1.2, alpha=0.95),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.2', color='black', lw=1.0)
            )

        # Use Log scale for Affinity (nM)
        ax.set_yscale('log')
        
        # Threshold annotations
        ax.axvline(0.5, color='green', linestyle=':', linewidth=1.2, label='SB %EL (0.5%)')
        ax.axvline(2.0, color='red', linestyle='--', linewidth=1.2, label='WB %EL (2.0%)')
        ax.axhline(500, color='black', linestyle='-.', linewidth=0.9, alpha=0.7, label='Affinity (500 nM)')
        
        ax.set_title(title_text, fontweight='bold')
        ax.set_xlabel('% EL Rank')
        ax.set_ylabel('Binding Affinity (nM) [Log Scale]')
        ax.legend(loc='lower right', frameon=True, fontsize=8)

    # Panel A: 1-AA (Affinity vs %EL Rank)
    ax1 = fig.add_subplot(2, 2, 1)
    plot_arm_panel(ax1, '1AA', 'A. 1-AA Selected Variants', CATEGORY_COLORS['1AA'])

    # Panel B: 2-AA (Affinity vs %EL Rank)
    ax2 = fig.add_subplot(2, 2, 2)
    plot_arm_panel(ax2, '2AA', 'B. 2-AA Selected Variants', CATEGORY_COLORS['2AA'])

    # Panel C: PS-CPL (Affinity vs %EL Rank)
    ax3 = fig.add_subplot(2, 2, 3)
    plot_arm_panel(ax3, 'PS-CPL', 'C. PS-CPL Selected Variants', CATEGORY_COLORS['PS-CPL'])

    # Panel D: 3D Comparative Distribution with Log10 Transformed CD8 Composite
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    
    for cat_name, color in CATEGORY_COLORS.items():
        sub_sel = selected_df[selected_df['category'].str.upper().str.replace('-', '') == cat_name.replace('-', '')].copy()
        if not sub_sel.empty:
            x = sub_sel[el_col]
            y = np.log10(np.clip(sub_sel[aff_col], 0.1, None))
            z = np.log10(np.clip(sub_sel['CD8_Composite'], 0, None) + 1.0)
            
            ax4.scatter(x, y, z, c=color, s=35, depthshade=True, edgecolor='k', linewidth=0.5, label=cat_name)

    ax4.set_title('D. Selected Cohort Profile (~150 Variants)', fontweight='bold', pad=10)
    ax4.set_xlabel('% EL Rank', labelpad=5)
    ax4.set_ylabel('Log10 Aff (nM)', labelpad=5)
    ax4.set_zlabel('Log10 (CD8 Comp + 1)', labelpad=5)
    ax4.legend(loc='upper left', frameon=True)
    ax4.view_init(elev=20, azim=-45)

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Dashboard saved successfully to {out_path}")
