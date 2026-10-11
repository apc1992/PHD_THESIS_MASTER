import json
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


def parse_pae_json_str(val) -> list:
  if isinstance(val, list):
    return val
  try:
    return json.loads(val)
  except Exception:
    return []


def main(
    pae_csv: Path = Path("output/mimotope_screening_pae_compiled.csv"),
    plddt_csv: Path = Path("output/peptide_anchor_stability.csv"),
    output_dir: Path = Path("output"),
    figure_dir: Path = Path("figures"),
    max_aggregator: float = 2.5,
    min_plddt: float = 80.0,
    max_anchor_pae_tier1: float = 4.0,
    max_anchor_pae_tier2: float = 4.8,
):
  if not pae_csv.exists() or not plddt_csv.exists():
    raise FileNotFoundError("Missing required input CSVs (PAE or pLDDT).")

  pae_df = pd.read_csv(pae_csv)
  plddt_df = pd.read_csv(plddt_csv)

  # 1. Summarize pLDDT per job/candidate
  group_col = "job_name" if "job_name" in plddt_df.columns else "peptide"
  plddt_grouped = []

  for key, group in plddt_df.groupby(group_col):
    mean_plddt = group["plddt"].mean()
    p2_row = group[group["residue_position"] == "P2"]
    p2_plddt = p2_row["plddt"].values[0] if not p2_row.empty else np.nan

    seq_len = group["length"].max() if "length" in group.columns else 9
    p9_row = group[group["residue_position"] == f"P{seq_len}"]
    p9_plddt = p9_row["plddt"].values[0] if not p9_row.empty else np.nan

    plddt_grouped.append({
        group_col: key,
        "mean_plddt": mean_plddt,
        "p2_plddt": p2_plddt,
        "p9_plddt": p9_plddt,
    })

  plddt_summary = pd.DataFrame(plddt_grouped)

  # 2. Merge PAE and pLDDT datasets
  merge_key = (
      "job_name"
      if "job_name" in pae_df.columns and "job_name" in plddt_summary.columns
      else "peptide"
  )
  merged_df = pd.merge(pae_df, plddt_summary, on=merge_key, how="inner")

  # Filter out controls from candidate selection
  if "category" in merged_df.columns:
    merged_df = merged_df.dropna(subset=["category"]).copy()
    merged_df["category_str"] = (
        merged_df["category"].astype(str).str.strip().str.lower()
    )
    control_terms = ["wt", "c57m", "neg", "control", "ctrl", "uncategorized"]
    merged_df = merged_df[~merged_df["category_str"].isin(control_terms)].copy()

  # 3. Extract anchor PAE scores
  merged_df["pae_list"] = merged_df["residue_pae_scores"].apply(
      parse_pae_json_str
  )
  merged_df["p2_pae"] = merged_df["pae_list"].apply(
      lambda x: x[1] if len(x) > 1 else np.nan
  )
  merged_df["p9_pae"] = merged_df["pae_list"].apply(
      lambda x: x[-1] if len(x) > 0 else np.nan
  )

  # 4. Multi-parametric selection tier logic
  pass_aggregator = merged_df["phla_aggregator_score"] <= max_aggregator
  pass_plddt = merged_df["mean_plddt"] >= min_plddt
  pass_tier1_anchors = (merged_df["p2_pae"] <= max_anchor_pae_tier1) & (
      merged_df["p9_pae"] <= max_anchor_pae_tier1
  )
  pass_tier2_anchors = (merged_df["p2_pae"] <= max_anchor_pae_tier2) & (
      merged_df["p9_pae"] <= max_anchor_pae_tier2
  )

  merged_df["rosetta_selection"] = "Filtered"
  merged_df.loc[
      pass_aggregator & pass_plddt & pass_tier2_anchors, "rosetta_selection"
  ] = "Tier 2"
  merged_df.loc[
      pass_aggregator & pass_plddt & pass_tier1_anchors, "rosetta_selection"
  ] = "Tier 1"

  # Export selected candidates
  selected_df = merged_df[
      merged_df["rosetta_selection"].isin(["Tier 1", "Tier 2"])
  ].copy()
  output_dir.mkdir(parents=True, exist_ok=True)
  out_csv = output_dir / "rosetta_candidates_selected.csv"
  selected_df.to_csv(out_csv, index=False)
  logging.info(
      f"Selected {len(selected_df)} candidates for Rosetta -> {out_csv}"
  )

  # 5. Publication Figure (2-Panel Layout)
  figure_dir.mkdir(parents=True, exist_ok=True)
  sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
  fig, axes = plt.subplots(
      1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [2, 1]}
  )

  # Panel A: Decision Matrix Scatter
  ax_a = axes[0]
  ax_a.axvspan(
      merged_df["phla_aggregator_score"].min() - 0.2,
      max_aggregator,
      min_plddt,
      100,
      color="#e6f2ff",
      alpha=0.5,
      label="Target Golden Zone",
  )
  ax_a.axvline(
      x=max_aggregator,
      color="#d95f02",
      linestyle="--",
      linewidth=1.5,
      alpha=0.8,
      label=f"PAE Aggregator Threshold (≤ {max_aggregator})",
  )
  ax_a.axhline(
      y=min_plddt,
      color="#7570b3",
      linestyle="--",
      linewidth=1.5,
      alpha=0.8,
      label=f"Peptide pLDDT Threshold (≥ {min_plddt})",
  )

  # Dynamic tier labels including anchor PAE thresholds
  tier_configs = {
      "Tier 1": (
          f"Selected Candidates - Tier 1 (anchors PAE ≤"
          f" {max_anchor_pae_tier1:.1f} Å)",
          "#1b9e77",
      ),
      "Tier 2": (
          f"Selected Candidates - Tier 2 (anchors PAE ≤"
          f" {max_anchor_pae_tier2:.1f} Å)",
          "#fdb462",
      ),
      "Filtered": (
          f"Filtered Out / Non-Binders (anchors PAE >"
          f" {max_anchor_pae_tier2:.1f} Å)",
          "#b0b0b0",
      ),
  }

  for status, (label_str, color) in tier_configs.items():
    sub = merged_df[merged_df["rosetta_selection"] == status]
    ax_a.scatter(
        sub["phla_aggregator_score"],
        sub["mean_plddt"],
        c=color,
        label=f"{label_str} (n={len(sub)})",
        alpha=0.85,
        s=70 if status != "Filtered" else 40,
        edgecolors="black" if status != "Filtered" else "none",
        linewidth=0.5,
    )

  ax_a.set_title(
      "A. Structural Selection Decision Matrix", fontsize=12, weight="bold"
  )
  ax_a.set_xlabel(
      "pHLA Aggregator Score (Lower = More Confident Interface Packing)",
      fontsize=10,
      weight="bold",
  )
  ax_a.set_ylabel(
      "Mean Peptide pLDDT (Higher = More Rigid Structural Backbone)",
      fontsize=10,
      weight="bold",
  )
  ax_a.legend(loc="upper right", frameon=True, facecolor="white", fontsize=8.5)

  # Panel B: Selection Yield by Design Category
  ax_b = axes[1]
  yield_df = (
      merged_df.groupby(["category", "rosetta_selection"])
      .size()
      .unstack(fill_value=0)
  )
  for col in ["Tier 1", "Tier 2", "Filtered"]:
    if col not in yield_df.columns:
      yield_df[col] = 0
  yield_df = yield_df[["Tier 1", "Tier 2", "Filtered"]]

  yield_df.plot(
      kind="bar",
      stacked=True,
      color=[
          tier_configs["Tier 1"][1],
          tier_configs["Tier 2"][1],
          tier_configs["Filtered"][1],
      ],
      ax=ax_b,
      edgecolor="black",
      linewidth=0.5,
  )
  ax_b.set_title(
      "B. Candidate Selection Yield by Category", fontsize=12, weight="bold"
  )
  ax_b.set_xlabel("Antigen Design Category", fontsize=10, weight="bold")
  ax_b.set_ylabel("Number of Candidates", fontsize=10, weight="bold")
  ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=0, ha="center")
  ax_b.legend(title="Status", frameon=True, facecolor="white", fontsize=9)

  plt.tight_layout()
  out_png = figure_dir / "rosetta_selection_decision_matrix.png"
  out_svg = figure_dir / "rosetta_selection_decision_matrix.svg"
  plt.savefig(out_png, dpi=300, bbox_inches="tight")
  plt.savefig(out_svg, bbox_inches="tight")
  plt.close()
  logging.info(
      f"Exported publication decision figures to {out_png} and {out_svg}"
  )


if __name__ == "__main__":
  main()
