from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def plot_category_heatmap(df: pd.DataFrame, category: str, output_dir: Path):
    sub_df = df[df["category"] == category].copy()
    n_peptides = sub_df["peptide"].nunique()
    
    if n_peptides == 0:
        print(f"[WARNING] No peptides found for category: {category}")
        return

    # Pivot data and compute mean pLDDT per residue position
    pivot_df = sub_df.pivot_table(
        index="residue_position",
        columns="peptide",
        values="plddt",
        aggfunc="mean"
    )

    # Order Y-axis residues explicitly P1..P9/P10
    y_order = [f"P{i}" for i in range(1, 11)]
    present_y = [p for p in y_order if p in pivot_df.index]
    pivot_df = pivot_df.reindex(present_y)

    # Remove axis index names to prevent seaborn labels on tick margins
    pivot_df.index.name = None
    pivot_df.columns.name = None

    # Dynamic figure width: scale with number of peptides
    fig_width = max(12, n_peptides * 0.28)
    fig, ax = plt.subplots(figsize=(fig_width, 5.5))

    # Heatmap with tuned contrast (vmin=70 to highlight position-specific variation)
    sns.heatmap(
        pivot_df,
        cmap="YlGnBu",
        vmin=70,
        vmax=100,
        cbar_kws={"label": "Residue pLDDT Score"},
        linewidths=0.5,
        linecolor="#ffffff",
        ax=ax,
    )

    ax.set_title(f"Position-Specific Residue Stability Landscape — Category: {category} (N={n_peptides})", fontsize=12, weight="bold", pad=12)
    ax.set_xlabel("Peptide Sequence", fontsize=10, weight="bold", labelpad=8)
    ax.set_ylabel("Residue Position", fontsize=10, weight="bold", labelpad=8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right", fontsize=9)

    out_file = output_dir / f"plddt_heatmap_{category.replace(' ', '_')}.png"
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SUCCESS] Exported {category} heatmap ({n_peptides} candidates) to: {out_file}")

def main(
    csv_path: Path = Path("output/peptide_anchor_stability.csv"),
    output_dir: Path = Path("figures"),
):
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # Filter out missing categories and standard controls
    df = df.dropna(subset=["category"]).copy()
    df["category_str"] = df["category"].astype(str).str.strip().str.lower()
    
    control_terms = ["wt", "c57m", "neg", "control", "ctrl", "uncategorized"]
    df = df[~df["category_str"].isin(control_terms)].copy()

    categories = sorted(df["category"].unique())
    output_dir.mkdir(parents=True, exist_ok=True)

    for cat in categories:
        plot_category_heatmap(df, cat, output_dir)

if __name__ == "__main__":
    main()
