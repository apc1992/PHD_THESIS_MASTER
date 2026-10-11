#!/usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd

# Set headless backend for HPC environments
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

WT_SEQUENCE = "QCFFCFKEL"

def parse_pae_string(pae_str):
    """Extracts float values from text string representation of PAE arrays."""
    if pd.isna(pae_str):
        return []
    vals = re.findall(r"[-+]?\d*\.\d+|\d+", str(pae_str))
    return [float(v) for v in vals]

def load_and_clean_master():
    master_path = "master_34_candidates_annotated.csv"
    if not os.path.exists(master_path):
        raise FileNotFoundError(f"Missing master file: {master_path}")
    
    df = pd.read_csv(master_path, low_memory=False)

    # Harmonize duplicate NetMHCpan columns
    if "HLA-A*02:01_%Rank_EL" in df.columns and "EL_rank_A0201" in df.columns:
        df["HLA-A*02:01_%Rank_EL"] = df["HLA-A*02:01_%Rank_EL"].fillna(df["EL_rank_A0201"])
    if "HLA-A*02:01_Aff(nM)" in df.columns and "BA_nM_A0201" in df.columns:
        df["HLA-A*02:01_Aff(nM)"] = df["HLA-A*02:01_Aff(nM)"].fillna(df["BA_nM_A0201"])

    # Drop non-target HLAs and duplicate legacy columns
    cols_to_drop = ["HLA-B*07:02_%Rank_EL", "HLA-C*04:01_%Rank_EL", "EL_rank_A0201", "BA_nM_A0201"]
    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True, errors="ignore")
    return df

def select_and_rank_positives(df):
    print("=== Processing Positive Candidates ===")
    
    # 1. Hard filters
    mask = (df["HLA-A*02:01_%Rank_EL"] < 2.0) & (df["HLA-A*02:01_Aff(nM)"] < 100.0)
    pos_df = df[mask].copy()

    # 2. Compute Z-Scores for Orthogonal Pillars (Lower Z = Better)
    z_rank = (pos_df["HLA-A*02:01_%Rank_EL"] - pos_df["HLA-A*02:01_%Rank_EL"].mean()) / (pos_df["HLA-A*02:01_%Rank_EL"].std() + 1e-6)
    
    pae_col = "raw_inter_pae" if "raw_inter_pae" in pos_df.columns else "phla_aggregator_score"
    z_af3 = (pos_df[pae_col] - pos_df[pae_col].mean()) / (pos_df[pae_col].std() + 1e-6)

    dg_col = "rosetta_dG_norm" if "rosetta_dG_norm" in pos_df.columns else "rosetta_dG_separated"
    dg = pos_df[dg_col].fillna(pos_df[dg_col].median())
    z_rosetta = (dg - dg.mean()) / (dg.std() + 1e-6)

    if "HLA-A*02:01_Thalf(h)" in pos_df.columns:
        thalf = pos_df["HLA-A*02:01_Thalf(h)"].fillna(pos_df["HLA-A*02:01_Thalf(h)"].median())
        z_thalf = (thalf - thalf.mean()) / (thalf.std() + 1e-6)
        pos_df["composite_rank_score"] = (0.30 * z_rank) + (0.30 * z_af3) + (0.25 * z_rosetta) - (0.15 * z_thalf)
    else:
        pos_df["composite_rank_score"] = (0.35 * z_rank) + (0.35 * z_af3) + (0.30 * z_rosetta)

    pos_df.sort_values(by="composite_rank_score", ascending=True, inplace=True)
    pos_df.reset_index(drop=True, inplace=True)

    print(f"[+] Selected {len(pos_df)} ranked positive candidates.")
    return pos_df

def generate_minimalist_plots(master_df, pos_ranked):
    """Generates publication-ready figures for Q1 journals and PPT slides."""
    print("=== Generating Minimalist Publication Plots ===")
    plot_dir = "final_selection_summary/plots"
    os.makedirs(plot_dir, exist_ok=True)

    # Global Style Parameters for Nature/Cell Minimalist Aesthetic
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.edgecolor": "#333333",
        "axes.linewidth": 1.0,
        "grid.color": "#e5e5e5",
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.color": "#333333",
        "ytick.color": "#333333"
    })

    # Prepare selection column
    pos_peptides = set(pos_ranked["peptide"] if "peptide" in pos_ranked.columns else pos_ranked["Peptide"])
    pep_col = "peptide" if "peptide" in master_df.columns else "Peptide"
    master_df["Selection_Status"] = master_df[pep_col].apply(
        lambda x: "Selected Positives" if x in pos_peptides else "Excluded / Filtered Out"
    )

    # -------------------------------------------------------------------------
    # FIGURE 1: Dual Quality Gate Scatter Plot (Affinity vs %Rank_EL)
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    rank_col = "HLA-A*02:01_%Rank_EL"
    aff_col = "HLA-A*02:01_Aff(nM)"

    sns.scatterplot(
        data=master_df,
        x=rank_col,
        y=aff_col,
        hue="Selection_Status",
        palette={"Selected Positives": "#1b9e77", "Excluded / Filtered Out": "#d95f02"},
        style="Selection_Status",
        markers={"Selected Positives": "o", "Excluded / Filtered Out": "X"},
        s=80,
        alpha=0.9,
        ax=ax
    )

    # Cutoff Threshold Lines
    ax.axvline(2.0, color="#666666", linestyle="--", linewidth=1.2, label="Rank Gate (2.0%)")
    ax.axhline(100.0, color="#666666", linestyle=":", linewidth=1.2, label="Affinity Gate (100 nM)")

    # Labeling and Formatting
    ax.set_xlabel("HLA-A*02:01 Elution Rank (%Rank EL)", fontsize=11, fontweight="bold")
    ax.set_ylabel("HLA-A*02:01 Binding Affinity (nM)", fontsize=11, fontweight="bold")
    ax.set_title("Tier 1: Positive Candidate Selection Quadrant", fontsize=12, fontweight="bold", pad=12)
    
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=True, facecolor="white", edgecolor="none", fontsize=9)

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "Figure1_Selection_Quadrant.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "Figure1_Selection_Quadrant.pdf"))
    plt.close()

    # -------------------------------------------------------------------------
    # FIGURE 2: Ranked Composite Score Bar Chart for Positives
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=300)

    # Identify peptide name column
    pos_pep_col = "peptide" if "peptide" in pos_ranked.columns else "Peptide"
    
    # Top 5 highlighted color palette
    colors = ["#08519c" if i < 5 else "#6baed6" for i in range(len(pos_ranked))]

    bars = ax.bar(
        range(len(pos_ranked)),
        pos_ranked["composite_rank_score"],
        color=colors,
        width=0.7,
        edgecolor="#333333",
        linewidth=0.5
    )

    ax.set_xticks(range(len(pos_ranked)))
    ax.set_xticklabels(pos_ranked[pos_pep_col], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Candidate Peptides (Ranked)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Composite Score (Lower = Better)", fontsize=11, fontweight="bold")
    ax.set_title("Tier 2: Orthogonal 4-Pillar Ranking of Positives", fontsize=12, fontweight="bold", pad=12)
    
    # Reference Line at zero
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="-")
    
    # Annotate Top 5
    ax.text(2, pos_ranked["composite_rank_score"].min() - 0.2, "Top 5 Candidates", 
            ha="center", fontsize=9, fontweight="bold", color="#08519c")

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "Figure2_Composite_Ranking.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "Figure2_Composite_Ranking.pdf"))
    plt.close()

    print(f"[+] Plots successfully saved to '{plot_dir}/'")

def select_unbiased_binding_negatives():
    print("=== Processing Binding/Cytokine Negative Controls ===")
    sources = [
        "selection_outputs/final_af3_2aa_candidates.csv",
        "selection_outputs/prioritized_leads_1AA.csv",
        "step1_all_variants.csv"
    ]
    
    dfs = []
    for src in sources:
        if os.path.exists(src):
            tmp = pd.read_csv(src, low_memory=False)
            tmp.rename(columns={"Peptide": "peptide", "EL_rank_A0201": "EL_rank", "BA_nM_A0201": "BA_nM"}, inplace=True, errors="ignore")
            dfs.append(tmp)
            
    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["peptide"])
    
    rank_col = "HLA-A*02:01_%Rank_EL" if "HLA-A*02:01_%Rank_EL" in combined.columns else "EL_rank"
    aff_col = "HLA-A*02:01_Aff(nM)" if "HLA-A*02:01_Aff(nM)" in combined.columns else "BA_nM"
    
    non_binders = combined[(combined[rank_col] > 10.0) | (combined[aff_col] > 10000.0)].copy()

    def count_subs(pep):
        if not isinstance(pep, str) or len(pep) != len(WT_SEQUENCE):
            return np.nan
        return sum(1 for a, b in zip(pep, WT_SEQUENCE) if a != b)

    if "num_substitutions" not in non_binders.columns or non_binders["num_substitutions"].isna().all():
        non_binders["num_substitutions"] = non_binders["peptide"].apply(count_subs)

    sampled_negs = []
    for subs, group in non_binders.groupby("num_substitutions"):
        sorted_g = group.sort_values(by=rank_col, ascending=False)
        sampled_negs.append(sorted_g.head(2))

    result = pd.concat(sampled_negs, ignore_index=True).sort_values(by=rank_col, ascending=False).head(5)
    
    keep_cols = ["peptide", "num_substitutions", rank_col, aff_col, "CD8_Composite", "Selection_Tier"]
    existing_keep = [c for c in keep_cols if c in result.columns]
    print(f"[+] Selected {len(result)} unbiased binding negative controls.")
    return result[existing_keep]

def select_structural_negatives():
    print("=== Processing Structural Negative Controls ===")
    pae_path = "af3/selection_outputs/mimotope_screening_best_models_compiled.csv"
    plddt_path = "af3/selection_outputs/peptide_anchor_stability.csv"
    
    if not (os.path.exists(pae_path) and os.path.exists(plddt_path)):
        print("[-] Warning: AF3 detailed structural files not found.")
        return pd.DataFrame()

    pae_df = pd.read_csv(pae_path)
    plddt_df = pd.read_csv(plddt_path)

    plddt_summary = []
    for pep, group in plddt_df.groupby("peptide"):
        mean_plddt = group["plddt"].mean()
        p2_row = group[group["residue_position"] == "P2"]
        p2_plddt = p2_row["plddt"].values[0] if not p2_row.empty else np.nan
        
        p9_row = group[group["anchor_type"] == "P-Omega_Anchor"]
        if p9_row.empty:
            max_len = group["length"].max()
            p9_row = group[group["residue_position"] == f"P{max_len}"]
        p9_plddt = p9_row["plddt"].values[0] if not p9_row.empty else np.nan
        
        plddt_summary.append({
            "peptide": pep,
            "mean_plddt": mean_plddt,
            "p2_plddt": p2_plddt,
            "p9_plddt": p9_plddt
        })
    plddt_sum_df = pd.DataFrame(plddt_summary)

    struct_df = pd.merge(pae_df, plddt_sum_df, on="peptide", how="inner")
    
    struct_df["pae_list"] = struct_df["residue_pae_scores"].apply(parse_pae_string)
    struct_df["p2_pae"] = struct_df["pae_list"].apply(lambda x: x[1] if len(x) > 1 else np.nan)
    struct_df["p9_pae"] = struct_df["pae_list"].apply(lambda x: x[-1] if len(x) > 0 else np.nan)
    struct_df.drop(columns=["pae_list", "residue_pae_scores"], inplace=True, errors="ignore")

    struct_df["worst_structural_score"] = struct_df["phla_aggregator_score"] + struct_df["p2_pae"] + struct_df["p9_pae"] - (0.1 * struct_df["mean_plddt"])
    worst_struct = struct_df.sort_values(by="worst_structural_score", ascending=False).head(5).copy()
    worst_struct.drop(columns=["worst_structural_score"], inplace=True)

    print(f"[+] Selected {len(worst_struct)} detailed structural negative controls.")
    return worst_struct

def main():
    os.makedirs("final_selection_summary", exist_ok=True)
    
    master_df = load_and_clean_master()
    pos_ranked = select_and_rank_positives(master_df)
    bind_negs = select_unbiased_binding_negatives()
    struct_negs = select_structural_negatives()

    # Generate Minimalist Plots
    generate_minimalist_plots(master_df, pos_ranked)

    # Save finalized outputs
    pos_ranked.to_csv("final_selection_summary/final_positive_candidates_ranked.csv", index=False)
    bind_negs.to_csv("final_selection_summary/negative_controls_binding_cytokine.csv", index=False)
    struct_negs.to_csv("final_selection_summary/negative_controls_structural.csv", index=False)

    print("\nSuccessfully updated files in 'final_selection_summary/':")
    print(f" - final_positive_candidates_ranked.csv ({len(pos_ranked)} candidates)")
    print(f" - negative_controls_binding_cytokine.csv ({len(bind_negs)} controls)")
    print(f" - negative_controls_structural.csv ({len(struct_negs)} controls)")
    print(f" - plots/Figure1_Selection_Quadrant.[png/pdf]")
    print(f" - plots/Figure2_Composite_Ranking.[png/pdf]")

if __name__ == "__main__":
    main()
