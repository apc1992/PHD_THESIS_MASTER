# pipeline/04_data_summary/scripts/publication_plots.py

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Standard 1Q Publication Pastel Palette across category metrics
CATEGORY_PASTEL = {
    '1AA': '#7eb0d5',
    '2AA': '#b2e061',
    'PS-CPL': '#ffb55a',
    'PS-CPL guided': '#fd7f6f'
}

# Unified Color Palette for Binder / Stability Status
STATUS_PALETTE = {
    'Strong Binder': '#d7191c', 'Strong Stable': '#d7191c',
    'Weak Binder': '#fdae61',   'Weak Stable': '#fdae61',
    'Non-Binder': '#d9d9d9',    'Non-Stable': '#d9d9d9'
}

def generate_binding_summary_figure(df: pd.DataFrame, output_dir: Path, peptide_name: str):
    """Generates 4-panel binding and stability diagnostics with enhanced typography and unified palettes."""
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.3)
    fig, axes = plt.subplots(2, 2, figsize=(16, 13), constrained_layout=True)
    axes = axes.flatten()

    # Filter out wild-type baseline
    df_var = df[df['category'] != 'WT'].copy()

    el_col = [c for c in df_var.columns if '%Rank_EL' in c][0]
    ba_col = [c for c in df_var.columns if 'Aff(nM)' in c][0]
    stab_rank_col = [c for c in df_var.columns if '%Rank_Stab' in c][0]
    thalf_col = [c for c in df_var.columns if 'Thalf(h)' in c][0]

    # Categorize Binders
    def classify_binder(row):
        el = row[el_col]
        aff = row[ba_col]
        if el <= 0.5 or aff <= 100:
            return 'Strong Binder'
        elif el <= 2.0 or aff <= 500:
            return 'Weak Binder'
        return 'Non-Binder'

    df_var['Binder_Class'] = df_var.apply(classify_binder, axis=1)

    # Categorize Stability (Fixed Thresholds: <= 0.5% SB, <= 2.0% WB)
    def classify_stability(row):
        stab = row[stab_rank_col]
        if stab <= 0.5:
            return 'Strong Stable'
        elif stab <= 2.0:
            return 'Weak Stable'
        return 'Non-Stable'

    df_var['Stability_Class'] = df_var.apply(classify_stability, axis=1)

    # Panel A: EL Rank vs. Affinity Scatter Plot
    sns.scatterplot(
        data=df_var, x=el_col, y=ba_col, hue='Binder_Class', 
        palette=STATUS_PALETTE, ax=axes[0], alpha=0.6, s=25
    )
    axes[0].set_yscale('log')
    axes[0].axvline(0.5, color='#d7191c', linestyle='--', lw=1.5, label='SB EL (<=0.5%)')
    axes[0].axvline(2.0, color='#fdae61', linestyle='--', lw=1.5, label='WB EL (<=2.0%)')
    axes[0].axhline(100, color='#d7191c', linestyle=':', lw=1.5, label='SB Aff (<=100nM)')
    axes[0].axhline(500, color='#fdae61', linestyle=':', lw=1.5, label='WB Aff (<=500nM)')
    axes[0].set_title("A. EL Rank vs. Binding Affinity (nM)", fontweight='bold', fontsize=14)
    axes[0].set_xlabel("% EL Rank", fontsize=12)
    axes[0].set_ylabel("Affinity (nM)", fontsize=12)
    axes[0].legend(loc='upper right', fontsize=11, title_fontsize=12, frameon=True)

    # Panel B: Binder Classification Yields
    categories = [c for c in ['1AA', '2AA', 'PS-CPL', 'PS-CPL guided'] if c in df_var['category'].unique()]
    yield_data = []

    for cat in categories:
        sub = df_var[df_var['category'] == cat]
        yield_data.extend([
            {'Category': cat, 'Metric': 'Total', 'Count': len(sub)},
            {'Category': cat, 'Metric': 'SB (EL <= 0.5)', 'Count': (sub[el_col] <= 0.5).sum()},
            {'Category': cat, 'Metric': 'SB (Aff <= 100)', 'Count': (sub[ba_col] <= 100).sum()},
            {'Category': cat, 'Metric': 'SB (Both)', 'Count': ((sub[el_col] <= 0.5) & (sub[ba_col] <= 100)).sum()},
            {'Category': cat, 'Metric': 'WB (EL <= 2.0)', 'Count': (sub[el_col] <= 2.0).sum()},
            {'Category': cat, 'Metric': 'WB (Aff <= 500)', 'Count': (sub[ba_col] <= 500).sum()},
            {'Category': cat, 'Metric': 'WB (Both)', 'Count': ((sub[el_col] <= 2.0) & (sub[ba_col] <= 500)).sum()},
            {'Category': cat, 'Metric': 'Non-Binders', 'Count': (sub['Binder_Class'] == 'Non-Binder').sum()}
        ])

    yield_df = pd.DataFrame(yield_data)
    sns.barplot(data=yield_df, x='Metric', y='Count', hue='Category', ax=axes[1], palette=CATEGORY_PASTEL)
    axes[1].set_yscale('log')
    axes[1].set_title("B. Binder Classification Yields (Counts)", fontweight='bold', fontsize=14)
    axes[1].set_ylabel("Number of Peptides (Log Scale)", fontsize=12)
    axes[1].tick_params(axis='x', rotation=35, labelsize=10)
    axes[1].legend(loc='upper right', fontsize=11, title="Category")

    # Panel C: %Rank Stability vs. Thalf (hours)
    sns.scatterplot(
        data=df_var, x=stab_rank_col, y=thalf_col, hue='Stability_Class',
        palette=STATUS_PALETTE, ax=axes[2], alpha=0.6, s=25
    )
    axes[2].axvline(0.5, color='#d7191c', linestyle='--', lw=1.5, label='Strong Stable (<=0.5%)')
    axes[2].axvline(2.0, color='#fdae61', linestyle='--', lw=1.5, label='Weak Stable (<=2.0%)')
    axes[2].set_title("C. %Rank Stability vs. Half-Life (T1/2 hours)", fontweight='bold', fontsize=14)
    axes[2].set_xlabel("%Rank Stability", fontsize=12)
    axes[2].set_ylabel("Half-Life T1/2 (hours)", fontsize=12)
    axes[2].legend(loc='upper right', fontsize=11, title_fontsize=12, frameon=True)

    # Panel D: Stability Classification Yields
    stab_yield_data = []
    for cat in categories:
        sub = df_var[df_var['category'] == cat]
        stab_yield_data.extend([
            {'Category': cat, 'Metric': 'Total', 'Count': len(sub)},
            {'Category': cat, 'Metric': 'Strong Stable (<=0.5%)', 'Count': (sub[stab_rank_col] <= 0.5).sum()},
            {'Category': cat, 'Metric': 'Weak Stable (<=2.0%)', 'Count': (sub[stab_rank_col] <= 2.0).sum()},
            {'Category': cat, 'Metric': 'Non-Stables', 'Count': (sub[stab_rank_col] > 2.0).sum()}
        ])

    stab_yield_df = pd.DataFrame(stab_yield_data)
    sns.barplot(data=stab_yield_df, x='Metric', y='Count', hue='Category', ax=axes[3], palette=CATEGORY_PASTEL)
    axes[3].set_yscale('log')
    axes[3].set_title("D. Stability Classification Yields (Counts)", fontweight='bold', fontsize=14)
    axes[3].set_ylabel("Number of Peptides (Log Scale)", fontsize=12)
    axes[3].tick_params(axis='x', rotation=20, labelsize=10)
    axes[3].legend(loc='upper right', fontsize=11, title="Category")

    fig.suptitle(f"[{peptide_name}] pMHC Binding and Stability Diagnostics", fontsize=16, fontweight='bold')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / f"{peptide_name}_binding_predictions_summary.png", dpi=300)
    plt.close()


def generate_functional_immunogenicity_figure(df: pd.DataFrame, output_dir: Path, peptide_name: str):
    """Generates structured functional immunogenicity figures with stratified category panels."""
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    
    df_var = df[df['category'] != 'WT'].copy()
    categories = [c for c in ['1AA', '2AA', 'PS-CPL', 'PS-CPL guided'] if c in df_var['category'].unique()]
    
    el_col = [c for c in df_var.columns if '%Rank_EL' in c][0]
    cyto_score_cols = [c for c in df_var.columns if c.endswith('_score')]

    # Create master figure layout using GridSpec
    fig = plt.figure(figsize=(18, 16), constrained_layout=True)
    gs = gridspec.GridSpec(3, 3, figure=fig, height_ratios=[1, 1.2, 1])

    # -------------------------------------------------------------------------
    # Panel A: Stratified Neo-Score Landscape (1x3 Sub-panels)
    # -------------------------------------------------------------------------
    for idx, cat in enumerate(categories[:3]):
        ax_a = fig.add_subplot(gs[0, idx])
        sub = df_var[df_var['category'] == cat]
        sns.scatterplot(
            data=sub, x='Foreignness', y='Dissimilarity', 
            color=CATEGORY_PASTEL.get(cat, '#7eb0d5'), ax=ax_a, alpha=0.7, s=30
        )
        ax_a.set_title(f"A{idx+1}. Neo-Score Landscape ({cat})", fontweight='bold', fontsize=12)
        ax_a.set_xlim(-0.05, 1.05)
        ax_a.set_ylim(-0.05, 1.05)

    # -------------------------------------------------------------------------
    # Panel B: Cytokine Activity Distributions (Full Width Middle Left/Center)
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[1, :2])
    if cyto_score_cols:
        melted_cyto = df_var.melt(
            id_vars=['category'], value_vars=cyto_score_cols, 
            var_name='Cytokine', value_name='Relative_Score'
        )
        # Clean cytokine names for plotting
        melted_cyto['Cytokine'] = melted_cyto['Cytokine'].str.replace('_score', '')
        sns.boxplot(
            data=melted_cyto, x='Cytokine', y='Relative_Score', hue='category', 
            hue_order=categories, ax=ax_b, palette=CATEGORY_PASTEL
        )
        ax_b.set_yscale('log')
        ax_b.set_title("B. Cytokine Activity Distributions across Categories", fontweight='bold', fontsize=14)
        ax_b.set_ylabel("Relative Score vs. WT (Log Scale)", fontsize=12)
        ax_b.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), frameon=True, fontsize=10, title="Category")

    # -------------------------------------------------------------------------
    # Panel D: Multi-Metric Correlation Matrix (Middle Right)
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 2])
    ba_col = [c for c in df_var.columns if 'Aff(nM)' in c][0]
    thalf_col = [c for c in df_var.columns if 'Thalf(h)' in c][0]
    stab_rank_col = [c for c in df_var.columns if '%Rank_Stab' in c][0]

    corr_candidates = [el_col, ba_col, thalf_col, stab_rank_col, 'Foreignness', 'Dissimilarity'] + cyto_score_cols
    valid_corr_cols = [c for c in corr_candidates if c in df_var.columns]

    corr = df_var[valid_corr_cols].corr()
    
    # Explicitly label CD4 / CD8 cytokine columns
    clean_labels = []
    for c in valid_corr_cols:
        if c.endswith('_score'):
            clean_labels.append(c.replace('_score', ''))
        elif '%Rank_EL' in c:
            clean_labels.append('EL Rank')
        elif 'Aff(nM)' in c:
            clean_labels.append('Affinity')
        elif 'Thalf(h)' in c:
            clean_labels.append('Thalf')
        elif '%Rank_Stab' in c:
            clean_labels.append('Stab Rank')
        else:
            clean_labels.append(c)

    sns.heatmap(
        corr, annot=True, fmt='.2f', cmap="vlag", vmin=-1, vmax=1, 
        ax=ax_d, xticklabels=clean_labels, yticklabels=clean_labels,
        cbar_kws={'shrink': 0.8}
    )
    ax_d.set_title("D. Multi-Metric Correlation (Pearson r)", fontweight='bold', fontsize=13)
    ax_d.tick_params(axis='x', rotation=45, labelsize=9)
    ax_d.tick_params(axis='y', labelsize=9)

    # -------------------------------------------------------------------------
    # Panel C: Stratified Cytokine vs %Rank EL Trade-off (1x3 Sub-panels)
    # -------------------------------------------------------------------------
    cd8_col = 'CD8_IFNg_score' if 'CD8_IFNg_score' in df_var.columns else cyto_score_cols[0]
    for idx, cat in enumerate(categories[:3]):
        ax_c = fig.add_subplot(gs[2, idx])
        sub = df_var[df_var['category'] == cat]
        sns.scatterplot(
            data=sub, x=el_col, y=cd8_col, 
            color=CATEGORY_PASTEL.get(cat, '#7eb0d5'), ax=ax_c, alpha=0.6, s=30
        )
        ax_c.set_yscale('log')
        ax_c.set_title(f"C{idx+1}. CD8 IFNg vs. %Rank EL ({cat})", fontweight='bold', fontsize=12)
        ax_c.set_xlabel("% EL Rank", fontsize=11)
        ax_c.set_ylabel("CD8 IFNg Score (Log Scale)" if idx == 0 else "", fontsize=11)

    fig.suptitle(f"[{peptide_name}] Immunogenicity & Functional Summary", fontsize=16, fontweight='bold')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / f"{peptide_name}_functional_immunogenicity_summary.png", dpi=300)
    plt.close()
