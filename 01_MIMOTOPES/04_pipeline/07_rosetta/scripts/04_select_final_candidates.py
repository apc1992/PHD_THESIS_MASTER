import logging
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


def main(
    rosetta_summary_csv: Path = Path("output/rosetta_candidate_summary.csv"),
    output_dir: Path = Path("output"),
    figure_dir: Path = Path("figures"),
    max_dG: float = -20.0,
    min_sc: float = 0.60,
    top_n_per_category: int = 5,
):
  if not rosetta_summary_csv.exists():
    raise FileNotFoundError(
        f"Rosetta summary CSV not found: {rosetta_summary_csv}"
    )

  df = pd.read_csv(rosetta_summary_csv)

  # Filter control entries
  df = df.dropna(subset=["category"]).copy()
  df["category_str"] = df["category"].astype(str).str.strip().str.lower()
  control_terms = ["wt", "c57m", "neg", "control", "ctrl", "uncategorized"]
  df = df[~df["category_str"].isin(control_terms)].copy()

  # Multi-parametric candidate evaluation logic
  pass_dG = df["mean_dG"] <= max_dG
  pass_sc = df["mean_sc_value"] >= min_sc

  df["final_selection_status"] = "Filtered"
  df.loc[pass_dG & pass_sc, "final_selection_status"] = "Lead Candidate"
  df.loc[pass_dG & ~pass_sc, "final_selection_status"] = "Secondary Candidate"

  # Rank candidates by binding energy within each category
  df = df.sort_values(by=["category", "mean_dG"], ascending=[True, True])

  # Select top N candidates per category
  selected_df = (
      df[df["final_selection_status"] != "Filtered"]
      .groupby("category")
      .head(top_n_per_category)
      .copy()
  )

  output_dir.mkdir(parents=True, exist_ok=True)
  out_csv = output_dir / "rosetta_final_candidates.csv"
  selected_df.to_csv(out_csv, index=False)
  logging.info(
      f"Selected {len(selected_df)} final candidates across design categories"
      f" -> {out_csv}"
  )

  # Publication-grade Final Ranking Matrix Figure
  figure_dir.mkdir(parents=True, exist_ok=True)
  sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
  fig, ax = plt.subplots(figsize=(10, 6))

  status_colors = {
      "Lead Candidate": "#1b9e77",
      "Secondary Candidate": "#fdb462",
      "Filtered": "#b0b0b0",
  }

  sns.scatterplot(
      data=df,
      x="mean_sc_value",
      y="mean_dG",
      hue="final_selection_status",
      style="category",
      palette=status_colors,
      s=90,
      ax=ax,
      edgecolor="black",
      linewidth=0.5,
  )

  ax.axhline(
      y=max_dG,
      color="#d95f02",
      linestyle="--",
      linewidth=1.2,
      label=f"$\Delta G$ Threshold (≤ {max_dG})",
  )
  ax.axvline(
      x=min_sc,
      color="#7570b3",
      linestyle="--",
      linewidth=1.2,
      label=f"Shape Comp. Threshold (≥ {min_sc})",
  )

  ax.set_title(
      "Final Mimotope Selection Matrix (Rosetta Structural Refinement)",
      fontsize=12,
      weight="bold",
  )
  ax.set_xlabel(
      "Shape Complementarity ($S_c$)", fontsize=10, weight="bold"
  )
  ax.set_ylabel(
      "Interface Binding Energy ($\Delta G_{\\text{separated}}$, REU)",
      fontsize=10,
      weight="bold",
  )
  ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, facecolor="white", fontsize=9)

  plt.tight_layout()
  out_png = figure_dir / "final_candidate_ranking_matrix.png"
  out_svg = figure_dir / "final_candidate_ranking_matrix.svg"
  plt.savefig(out_png, dpi=300, bbox_inches="tight")
  plt.savefig(out_svg, bbox_inches="tight")
  plt.close()
  logging.info(f"Exported candidate ranking figures to {out_png} and {out_svg}")


if __name__ == "__main__":
  main()
