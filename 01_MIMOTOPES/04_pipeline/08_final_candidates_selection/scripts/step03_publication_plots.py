#!/usr/bin/env python3
"""step03_publication_plots.py

Generates publication-ready figures for Tier 1 quadrant selection and Tier 2 composite Z-score ranking.
"""

import logging
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger("08_final_candidates")


def plot_selection_funnel(df: pd.DataFrame, output_dir: Path):
  """Figure 1: Positive Candidate Selection Quadrant."""
  plt.figure(figsize=(8, 6), dpi=300)

  if df.empty:
    logger.warning("Master table is empty. Unable to plot selection funnel.")
    return

  passed = df[df["passed_all_gates"] == True]
  excluded = df[df["passed_all_gates"] == False]

  plt.scatter(
      excluded["elution_rank"],
      excluded["binding_affinity_nM"],
      color="#d95f02",
      marker="x",
      s=65,
      label=f"Excluded ({len(excluded)})",
      alpha=0.85,
      zorder=2,
  )
  plt.scatter(
      passed["elution_rank"],
      passed["binding_affinity_nM"],
      color="#1b9e77",
      marker="o",
      s=75,
      label=f"Selected Positives ({len(passed)})",
      alpha=0.9,
      zorder=3,
  )

  plt.yscale("log")
  plt.ylim(5, 25000)

  # Gates: 2.5% EL Rank and 500 nM Affinity
  plt.axvline(
      x=2.5,
      color="gray",
      linestyle="--",
      linewidth=1.2,
      label="Rank Gate (2.5%)",
      zorder=1,
  )
  plt.axhline(
      y=500.0,
      color="gray",
      linestyle=":",
      linewidth=1.2,
      label="Affinity Gate (500 nM)",
      zorder=1,
  )

  plt.xlabel(
      "HLA-A*02:01 Elution Rank (%Rank EL)", fontweight="bold", fontsize=11
  )
  plt.ylabel(
      "HLA-A*02:01 Binding Affinity (nM)", fontweight="bold", fontsize=11
  )
  plt.title(
      f"Tier 1: Positive Selection Quadrant (N={len(df)} AF3 Candidates)",
      fontweight="bold",
      fontsize=13,
      pad=12,
  )
  plt.legend(frameon=True, loc="lower right", fontsize=10)

  sns.despine(top=True, right=True)
  plt.tight_layout()

  out_path = output_dir / "candidate_selection_funnel.png"
  plt.savefig(out_path, bbox_inches="tight")
  plt.close()
  logger.info(f"Saved Figure 1 to {out_path}")


def plot_zscore_breakdown(
    df_selected_leads: pd.DataFrame, output_dir: Path, top_n: int = None
):
  """Figure 2: Orthogonal Composite Score Bar Plot of All Step 1 Passed Positives.

  Colors negative scores (better than average) in dark blue and positive scores
  in light blue.
  """
  df_plot = df_selected_leads.copy()
  if top_n is not None and len(df_plot) > top_n:
    df_plot = df_plot.head(top_n)

  if df_plot.empty:
    logger.warning("No candidates available to display in Figure 2.")
    return

  fig, ax = plt.subplots(figsize=(max(8, len(df_plot) * 0.4), 5.5), dpi=300)

  # Dark blue for negative composite scores (better), light blue for positive
  colors = [
      "#08519c" if score < 0 else "#6baed6" for score in df_plot["composite_score"]
  ]

  ax.bar(
      df_plot["peptide"],
      df_plot["composite_score"],
      color=colors,
      edgecolor="black",
      linewidth=0.6,
      width=0.7,
  )

  ax.axhline(0, color="black", linewidth=1.0, linestyle="-")
  ax.set_ylabel(
      "Composite Score (Lower = Better)", fontweight="bold", fontsize=11
  )
  ax.set_xlabel(
      "Selected Step 1 Candidates (Ranked)", fontweight="bold", fontsize=11
  )
  ax.set_title(
      f"Tier 2: Orthogonal 4-Pillar Ranking ({len(df_plot)} Selected Leads)",
      fontweight="bold",
      fontsize=13,
      pad=12,
  )

  plt.xticks(rotation=45, ha="right", fontsize=9, fontweight="bold")
  sns.despine(top=True, right=True)
  plt.tight_layout()

  out_path = output_dir / "lead_candidates_zscore_breakdown.png"
  plt.savefig(out_path, bbox_inches="tight")
  plt.close()
  logger.info(f"Saved Figure 2 to {out_path}")
