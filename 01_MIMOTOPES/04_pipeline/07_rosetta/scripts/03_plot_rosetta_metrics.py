import logging
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

try:
    from adjustText import adjust_text

    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


def main(
    summary_csv: Path = Path("output/rosetta_candidate_summary.csv"),
    compiled_csv: Path = Path("output/rosetta_compiled_scores.csv"),
    figure_dir: Path = Path("figures"),
):
    if not summary_csv.exists():
        raise FileNotFoundError(f"Missing summary CSV: {summary_csv}")

    summary_df = pd.read_csv(summary_csv)
    compiled_df = (
        pd.read_csv(compiled_csv) if compiled_csv.exists() else pd.DataFrame()
    )

    # Filter control artifacts
    summary_df = summary_df.dropna(subset=["category"]).copy()
    summary_df["category_str"] = (
        summary_df["category"].astype(str).str.strip().str.lower()
    )
    control_terms = ["wt", "c57m", "neg", "control", "ctrl", "uncategorized"]
    summary_df = summary_df[
        ~summary_df["category_str"].isin(control_terms)
    ].copy()

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    categories = sorted(summary_df["category"].unique())
    palette = sns.color_palette("Set2", n_colors=len(categories))
    cat_color_map = dict(zip(categories, palette))

    # Panel A: Interface Binding Energy (dG_separated) by Category
    sns.boxplot(
        data=summary_df,
        x="category",
        y="mean_dG",
        hue="category",
        order=categories,
        palette=palette,
        width=0.4,
        legend=False,
        ax=axes[0],
    )
    sns.stripplot(
        data=summary_df,
        x="category",
        y="mean_dG",
        order=categories,
        color="#2b2b2b",
        alpha=0.6,
        jitter=0.2,
        size=6,
        ax=axes[0],
    )
    axes[0].set_title(
        r"A. Interface Binding Energy ($\Delta G_{\text{separated}}$)",
        fontsize=11,
        weight="bold",
    )
    axes[0].set_ylabel(
        r"Mean $\Delta G_{\text{separated}}$ (REU)", fontsize=10, weight="bold"
    )
    axes[0].set_xlabel("Antigen Design Category", fontsize=10, weight="bold")
    axes[0].axhline(
        y=-25.0,
        color="#d95f02",
        linestyle="--",
        alpha=0.8,
        label=r"High Affinity Threshold ($\leq -25.0$)",
    )
    axes[0].legend(
        loc="lower left", frameon=True, facecolor="white", fontsize=8.5
    )

    # Panel B: Shape Complementarity (Sc) vs Normalized dG
    for cat in categories:
        sub = summary_df[summary_df["category"] == cat]
        axes[1].scatter(
            sub["mean_sc_value"],
            sub["mean_norm_dG"],
            color=cat_color_map[cat],
            label=cat,
            s=75,
            alpha=0.85,
            edgecolors="black",
            linewidth=0.5,
        )

    axes[1].axvline(
        x=0.65,
        color="#7570b3",
        linestyle="--",
        linewidth=1.2,
        label=r"Shape Comp. $S_c \geq 0.65$",
    )
    axes[1].set_title(
        r"B. Shape Complementarity ($S_c$) vs. Normalized $\Delta G$",
        fontsize=11,
        weight="bold",
    )
    axes[1].set_xlabel(
        r"Shape Complementarity ($S_c$)", fontsize=10, weight="bold"
    )
    axes[1].set_ylabel(
        r"$\Delta G_{\text{separated}} / \Delta SASA \times 100$",
        fontsize=10,
        weight="bold",
    )
    axes[1].legend(
        loc="upper left", frameon=True, facecolor="white", fontsize=8.5
    )

    # Annotate top-tier candidates (S_c >= 0.73) with repulsion lines
    top_candidates = summary_df[summary_df["mean_sc_value"] >= 0.73].copy()

    texts = []
    for _, row in top_candidates.iterrows():
        texts.append(
            axes[1].text(
                row["mean_sc_value"],
                row["mean_norm_dG"],
                row["peptide"],
                fontsize=7.0,
                alpha=0.9,
                weight="semibold",
            )
        )

    if HAS_ADJUST_TEXT and texts:
        adjust_text(
            texts,
            ax=axes[1],
            arrowprops=dict(
                arrowstyle="-", color="#444444", lw=0.5, alpha=0.6
            ),
            expand_points=(1.8, 1.8),
            expand_text=(1.4, 1.4),
            force_points=(0.3, 0.6),
            force_text=(0.5, 0.8),
        )

    # Panel C: Decoy Relaxation Energy Convergence across nstruct
    if not compiled_df.empty and "category" in compiled_df.columns:
        compiled_df = compiled_df[
            compiled_df["peptide"].isin(summary_df["peptide"])
        ].copy()
        sns.boxplot(
            data=compiled_df,
            x="category",
            y="total_score",
            hue="category",
            order=categories,
            palette=palette,
            width=0.4,
            legend=False,
            ax=axes[2],
        )
        sns.stripplot(
            data=compiled_df,
            x="category",
            y="total_score",
            order=categories,
            color="#2b2b2b",
            alpha=0.4,
            jitter=0.2,
            size=4,
            ax=axes[2],
        )
        axes[2].set_title(
            "C. Total Rosetta Energy Convergence Across Decoys",
            fontsize=11,
            weight="bold",
        )
        axes[2].set_ylabel(
            "Rosetta Total Score (REU)", fontsize=10, weight="bold"
        )
        axes[2].set_xlabel(
            "Antigen Design Category", fontsize=10, weight="bold"
        )

    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    out_png = figure_dir / "rosetta_interface_evaluation.png"
    out_svg = figure_dir / "rosetta_interface_evaluation.svg"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()
    logging.info(f"Saved Rosetta interface plots to {out_png} and {out_svg}")


if __name__ == "__main__":
    main()
