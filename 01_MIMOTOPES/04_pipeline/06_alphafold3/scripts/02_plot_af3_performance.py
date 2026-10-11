from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def main(
    input_csv: Path = Path("output/af3_compiled_results.csv"),
    output_png: Path = Path("figures/af3_category_performance.png"),
):
    if not input_csv.exists():
        raise FileNotFoundError(f"Compiled input CSV not found: {input_csv}")

    df = pd.read_csv(input_csv)

    # Clean and filter control entries
    df = df.dropna(subset=["category"]).copy()
    df["category_str"] = df["category"].astype(str).str.strip().str.lower()
    control_terms = ["wt", "c57m", "neg", "control", "ctrl", "uncategorized"]
    df = df[~df["category_str"].isin(control_terms)].copy()

    if df.empty:
        print("[WARNING] No non-control candidate data available to plot.")
        return

    # Style configuration
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    categories = sorted(df["category"].unique())
    palette = sns.color_palette("Set2", n_colors=len(categories))

    # Panel 1: ipTM Interface Distribution
    sns.violinplot(
        data=df, x="category", y="iptm", order=categories,
        palette=palette, inner="box", cut=0, ax=axes[0]
    )
    sns.stripplot(
        data=df, x="category", y="iptm", order=categories,
        color="#2b2b2b", alpha=0.5, jitter=0.2, size=5, ax=axes[0]
    )
    axes[0].set_title("Interface Confidence ($ipTM$) by Category", fontsize=12, weight="bold")
    axes[0].set_ylabel("AlphaFold 3 $ipTM$ Score", fontsize=10, weight="bold")
    axes[0].set_xlabel("Antigen Design Category", fontsize=10, weight="bold")
    axes[0].axhline(y=0.85, color="#d95f02", linestyle="--", alpha=0.8, label="High Quality (>0.85)")
    axes[0].axhline(y=0.70, color="#7570b3", linestyle="--", alpha=0.8, label="Acceptable (>0.70)")
    axes[0].legend(loc="lower left", frameon=True, facecolor="white", fontsize=9)

    # Panel 2: Total Ranking Score Distribution
    sns.boxplot(
        data=df, x="category", y="ranking_score", order=categories,
        palette=palette, width=0.4, ax=axes[1]
    )
    sns.stripplot(
        data=df, x="category", y="ranking_score", order=categories,
        color="#2b2b2b", alpha=0.5, jitter=0.2, size=5, ax=axes[1]
    )
    axes[1].set_title("Composite Ranking Score by Category", fontsize=12, weight="bold")
    axes[1].set_ylabel("Ranking Score", fontsize=10, weight="bold")
    axes[1].set_xlabel("Antigen Design Category", fontsize=10, weight="bold")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SUCCESS] Performance plot saved to: {output_png}")

if __name__ == "__main__":
    main()
