import os
import re
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ==========================================
# 1. Helper Functions & Normalization
# ==========================================


def extract_4digit_keys(locus_label, cell_value):
  """Converts locus row labels (e.g., 'A_1', 'DPA1_1', 'DRB1_2') and cell values

  (e.g., '02:01', 'DPA1*01:03:01') into normalized 4-digit keys (e.g., 'A0201',
  'DPA10103', 'DRB10101').
  """
  if pd.isna(cell_value) or not str(cell_value).strip():
    return set()

  gene = re.sub(r"_\d+$", "", str(locus_label).strip().upper())
  keys = set()

  for part in re.split(r"[,;]+", str(cell_value)):
    part = part.strip()
    if not part or part.lower() in ["nan", "none"]:
      continue

    part = re.sub(r"^HLA-", "", part, flags=re.IGNORECASE)
    if not part.upper().startswith(gene):
      part = f"{gene}*{part}"

    m = re.match(r"^([A-Z0-9]+)\*?(\d{2}):?(\d{2})", part, re.IGNORECASE)
    if m:
      g, f1, f2 = m.groups()
      keys.add(f"{g.upper()}{f1}{f2}")
    else:
      clean = (
          part.replace("*", "")
          .replace(":", "")
          .replace("_", "")
          .replace("-", "")
          .upper()
      )
      keys.add(clean)

  return keys


def get_required_allele_keys(mhc_raw):
  """Parses raw Class I or Class II MHC strings into required 4-digit allele keys."""
  mhc_str = str(mhc_raw).strip()

  if "-" in mhc_str or "/" in mhc_str:
    parts = re.split(r"[-/]+", mhc_str)
    keys = set()
    for p in parts:
      p_clean = (
          re.sub(r"^HLA-?", "", p, flags=re.IGNORECASE)
          .replace("*", "")
          .replace(":", "")
          .replace("_", "")
      )
      if p_clean:
        keys.add(p_clean)
    return keys
  else:
    p_clean = (
        re.sub(r"^HLA-?", "", mhc_str, flags=re.IGNORECASE)
        .replace("*", "")
        .replace(":", "")
        .replace("_", "")
    )
    return {p_clean}


def sort_key_hla_i(allele):
  """Sort key for Class I: A first, then B, then C."""
  m = re.search(r"HLA-([ABC])\*?(\d+):?(\d+)?", str(allele))
  if m:
    locus, g1, g2 = m.groups()
    locus_rank = {"A": 1, "B": 2, "C": 3}.get(locus, 9)
    return (locus_rank, int(g1), int(g2) if g2 else 0)
  return (99, 0, 0)


def sort_key_hla_ii(allele):
  """Sort key for Class II: DP first, then DQ, then DR."""
  s = str(allele).upper()
  if "DP" in s:
    locus_rank = 1
  elif "DQ" in s:
    locus_rank = 2
  elif "DR" in s:
    locus_rank = 3
  else:
    locus_rank = 4
  digits = [int(x) for x in re.findall(r"\d+", s)]
  return (locus_rank, digits)


def build_custom_cmap():
  """Creates custom colormap:

  - 0 to 5: darkred -> lightsalmon
  - 5 to 100: lightblue -> darkblue
  - Missing/Not Present: white
  """
  colors = [
      (0.0, "#8B0000"),  # Dark red at 0
      (0.05, "#FFA07A"),  # Light salmon at 5
      (0.05001, "#ADD8E6"),  # Light blue just above 5
      (1.0, "#00008B"),  # Dark blue at 100
  ]
  cmap = mcolors.LinearSegmentedColormap.from_list(
      "hla_rank_cmap", colors, N=256
  )
  cmap.set_bad("white")
  return cmap


# ==========================================
# 2. Typing Parser
# ==========================================


def load_leukopak_typing(typing_filepath):
  """Parses leukopaks_hla_typing.csv into a mapping:

  { leukopak_id: set_of_4digit_allele_keys }
  """
  df_typing = pd.read_csv(typing_filepath)

  locus_col = df_typing.columns[0]
  sample_cols = [c.strip() for c in df_typing.columns[1:]]

  typing_dict = {sample: set() for sample in sample_cols}

  for _, row in df_typing.iterrows():
    locus_label = str(row[locus_col]).strip()
    for sample in sample_cols:
      cell_val = row[sample]
      keys = extract_4digit_keys(locus_label, cell_val)
      typing_dict[sample].update(keys)

  return typing_dict


# ==========================================
# 3. Enhanced Multi-Panel Heatmap Plotter
# ==========================================


def plot_heatmap_with_counts(
    matrix, title, out_filename, cmap, rank_threshold=15.0
):
  """Renders heatmap alongside an aligned right-hand bar chart counting alleles

  with %Rank_EL < threshold per Leukopak.
  """
  n_rows, n_cols = matrix.shape

  # Dynamically scale figure dimensions for improved readability
  fig_width = max(16, n_cols * 0.55)
  fig_height = max(10, n_rows * 0.55)

  fig = plt.figure(figsize=(fig_width, fig_height))

  # Grid layout: [Heatmap, Colorbar, Side Bar Chart]
  gs = gridspec.GridSpec(1, 3, width_ratios=[14, 0.4, 2.5], wspace=0.08)

  ax_main = fig.add_subplot(gs[0])
  ax_cbar = fig.add_subplot(gs[1])
  ax_bar = fig.add_subplot(gs[2], sharey=ax_main)

  # --- 1. Main Heatmap ---
  sns.heatmap(
      matrix,
      ax=ax_main,
      cbar_ax=ax_cbar,
      cmap=cmap,
      vmin=0,
      vmax=100,
      annot=True,
      fmt=".1f",
      annot_kws={"size": 9.5, "weight": "bold"},
      linewidths=0.6,
      linecolor="#D3D3D3",
      cbar_kws={"label": "%Rank_EL / Rank_EL"},
  )

  ax_main.set_title(title, fontsize=16, pad=20, weight="bold")
  ax_main.set_xlabel("HLA Alleles", fontsize=13, labelpad=12, weight="bold")
  ax_main.set_ylabel("Leukopaks", fontsize=13, weight="bold")

  # Precise tick alignment anchor for rotated labels
  ax_main.tick_params(axis="x", labelsize=10.5)
  ax_main.tick_params(axis="y", labelsize=11, rotation=0)
  plt.setp(
      ax_main.get_xticklabels(),
      rotation=45,
      ha="right",
      rotation_mode="anchor",
  )

  # Colorbar typography
  ax_cbar.tick_params(labelsize=10)
  ax_cbar.yaxis.label.set_size(12)
  ax_cbar.yaxis.label.set_weight("bold")

  # --- 2. Right Side Bar Chart (< 15 count) ---
  counts = (matrix < rank_threshold).sum(axis=1).values
  y_positions = np.arange(n_rows) + 0.5

  bars = ax_bar.barh(
      y_positions,
      counts,
      height=0.65,
      color="#2b5c8f",
      edgecolor="#1a3654",
      align="center",
  )

  ax_bar.set_title(
      f"Alleles < {int(rank_threshold)}", fontsize=13, pad=20, weight="bold"
  )
  ax_bar.set_xlabel("Count", fontsize=12, weight="bold")
  ax_bar.tick_params(axis="x", labelsize=10)
  ax_bar.grid(axis="x", linestyle="--", alpha=0.5)

  # X-axis limits and value labels on bars
  max_c = max(counts) if len(counts) > 0 and max(counts) > 0 else 1
  ax_bar.set_xlim(0, max_c * 1.3 + 1)

  for bar, count in zip(bars, counts):
    if count > 0:
      ax_bar.text(
          bar.get_width() + 0.15,
          bar.get_y() + bar.get_height() / 2.0,
          f"{count}",
          va="center",
          ha="left",
          fontsize=10.5,
          weight="bold",
          color="#1a3654",
      )

  # Hide Y-axis labels on the bar plot panel
  plt.setp(ax_bar.get_yticklabels(), visible=False)

  plt.savefig(out_filename, dpi=300, bbox_inches="tight")
  plt.close()
  print(f"Saved: {out_filename}")


# ==========================================
# 4. Processing Loops
# ==========================================


def process_hla_class_i(class1_filepath, typing_dict, cmap):
  df_c1 = pd.read_csv(class1_filepath)

  rank_cols = [c for c in df_c1.columns if c.endswith("_%Rank_EL")]
  raw_alleles = [c.replace("_%Rank_EL", "") for c in rank_cols]
  sorted_alleles = sorted(raw_alleles, key=sort_key_hla_i)

  peptides = df_c1["Peptide"].unique()

  for pep in peptides:
    pep_df = df_c1[df_c1["Peptide"] == pep].iloc[0]
    matrix = pd.DataFrame(
        np.nan, index=list(typing_dict.keys()), columns=sorted_alleles
    )

    for sample_id, sample_keys in typing_dict.items():
      for raw_allele in sorted_alleles:
        req_keys = get_required_allele_keys(raw_allele)
        if req_keys.issubset(sample_keys):
          matrix.loc[sample_id, raw_allele] = pep_df[f"{raw_allele}_%Rank_EL"]

    plot_heatmap_with_counts(
        matrix=matrix,
        title=f"HLA Class I %Rank_EL - Peptide: {pep}",
        out_filename=f"HLA_I_heatmap_{pep}.png",
        cmap=cmap,
        rank_threshold=15.0,
    )


def process_hla_class_ii(dp_file, dq_file, drb_file, typing_dict, cmap):
  dfs = []
  for f in [dp_file, dq_file, drb_file]:
    if os.path.exists(f):
      dfs.append(pd.read_csv(f))

  if not dfs:
    print("No Class II prediction files found.")
    return

  df_c2 = pd.concat(dfs, ignore_index=True)
  pep_col = "mhc2_pep" if "mhc2_pep" in df_c2.columns else "Peptide"
  peptides = df_c2[pep_col].unique()

  unique_alleles = df_c2["MHC"].unique()
  sorted_alleles = sorted(unique_alleles, key=sort_key_hla_ii)

  for pep in peptides:
    pep_sub = df_c2[df_c2[pep_col] == pep]
    pred_lookup = dict(zip(pep_sub["MHC"], pep_sub["Rank_EL"]))

    matrix = pd.DataFrame(
        np.nan, index=list(typing_dict.keys()), columns=sorted_alleles
    )

    for sample_id, sample_keys in typing_dict.items():
      for raw_mhc in sorted_alleles:
        req_keys = get_required_allele_keys(raw_mhc)
        if req_keys.issubset(sample_keys) and raw_mhc in pred_lookup:
          matrix.loc[sample_id, raw_mhc] = pred_lookup[raw_mhc]

    plot_heatmap_with_counts(
        matrix=matrix,
        title=f"HLA Class II Rank_EL - Peptide: {pep}",
        out_filename=f"HLA_II_heatmap_{pep}.png",
        cmap=cmap,
        rank_threshold=15.0,
    )


# ==========================================
# Main Execution
# ==========================================

if __name__ == "__main__":
  TYPING_FILE = "leukopaks_hla_typing.csv"
  CLASS1_FILE = "netmhcpan_predictions_full_table.csv"
  DP_FILE = "dp_preds.csv"
  DQ_FILE = "dq_preds.csv"
  DRB_FILE = "drb_preds.csv"

  custom_cmap = build_custom_cmap()
  typing_data = load_leukopak_typing(TYPING_FILE)

  process_hla_class_i(CLASS1_FILE, typing_data, custom_cmap)
  process_hla_class_ii(DP_FILE, DQ_FILE, DRB_FILE, typing_data, custom_cmap)
