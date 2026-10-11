#!/usr/bin/env python3
import glob
import os
import pandas as pd

WT_PEPTIDE = "QCFFCFKEL"


def annotate_mutations(peptide, wt=WT_PEPTIDE):
    """Calculates number of substitutions, positions, and specific mutation strings."""
    if len(peptide) != len(wt):
        return len(peptide), "Length_Mismatch", "Length_Mismatch"

    mutations = []
    positions = []
    for i, (wt_aa, mut_aa) in enumerate(zip(wt, peptide)):
        if wt_aa != mut_aa:
            positions.append(i + 1)  # 1-based index
            mutations.append(f"{wt_aa}{i+1}{mut_aa}")

    num_subs = len(mutations)
    pos_str = ",".join(map(str, positions)) if positions else "None"
    mut_str = ", ".join(mutations) if mutations else "WT"
    return num_subs, pos_str, mut_str


def load_and_standardize(file_path):
    """Reads CSV with low_memory=False and renames 'Peptide' to 'peptide'."""
    if not os.path.exists(file_path):
        print(f"Warning: File not found -> {file_path}")
        return pd.DataFrame()
    df = pd.read_csv(file_path, low_memory=False)
    if "Peptide" in df.columns and "peptide" not in df.columns:
        df.rename(columns={"Peptide": "peptide"}, inplace=True)
    return df


def parse_rosetta_interface(rosetta_dir="rosetta/rosetta_outputs"):
    """Aggregates Rosetta InterfaceAnalyzer outputs across replicates."""
    score_files = glob.glob(os.path.join(rosetta_dir, "*/interface_score.sc"))
    rows = []
    for sf in score_files:
        pep = os.path.basename(os.path.dirname(sf))
        with open(sf, "r") as f:
            lines = f.readlines()

        header, data = None, []
        for line in lines:
            if line.startswith("SCORE:") and "total_score" in line:
                header = line.strip().split()[1:]
            elif line.startswith("SCORE:") and "total_score" not in line:
                data.append(line.strip().split()[1:])

        if header and data:
            for dl in data:
                if len(dl) == len(header):
                    r = dict(zip(header, dl))
                    r["peptide"] = pep
                    rows.append(r)

    if not rows:
        return pd.DataFrame()

    df_int = pd.DataFrame(rows)
    for col in [
        "dG_separated",
        "dG_separated/dSASAx100",
        "sc_value",
        "dSASA_int",
        "per_residue_energy_int",
    ]:
        if col in df_int.columns:
            df_int[col] = pd.to_numeric(df_int[col], errors="coerce")

    summary = (
        df_int.groupby("peptide")
        .agg(
            rosetta_dG_separated=("dG_separated", "mean"),
            rosetta_dG_norm=("dG_separated/dSASAx100", "mean"),
            rosetta_sc_value=("sc_value", "mean"),
            rosetta_dSASA_int=("dSASA_int", "mean"),
            rosetta_per_res_energy=("per_residue_energy_int", "mean"),
        )
        .reset_index()
    )
    return summary


def build_master_table():
    print("=== Building Master Candidate Table ===")

    # 1. Base 34 Candidate Peptides
    rosetta_sel_path = "af3/selection_outputs/rosetta_candidates_selected.csv"
    if os.path.exists(rosetta_sel_path):
        base_df = load_and_standardize(rosetta_sel_path)
    else:
        peps = [
            d
            for d in os.listdir("rosetta/rosetta_outputs")
            if os.path.isdir(os.path.join("rosetta/rosetta_outputs", d))
        ]
        base_df = pd.DataFrame({"peptide": peps})

    print(f"Loaded {len(base_df)} base candidates.")

    # 2. Annotate Mutations & Origin
    mut_info = base_df["peptide"].apply(lambda p: annotate_mutations(p))
    base_df["num_substitutions"] = [m[0] for m in mut_info]
    base_df["mut_positions"] = [m[1] for m in mut_info]
    base_df["mutations_vs_wt"] = [m[2] for m in mut_info]

    step2_leads = []
    if os.path.exists("selection_outputs/prioritized_leads_step2.csv"):
        step2_df = load_and_standardize("selection_outputs/prioritized_leads_step2.csv")
        if "peptide" in step2_df.columns:
            step2_leads = step2_df["peptide"].tolist()

    def assign_category(row):
        if row["peptide"] in step2_leads:
            return "Targeted (Step 2)"
        elif row["num_substitutions"] == 1:
            return "1AA"
        elif row["num_substitutions"] == 2:
            return "2AA"
        return "Other/Targeted"

    base_df["origin_source"] = base_df.apply(assign_category, axis=1)

    # 3. Merge Phase 1 Summary & Cytokine Metrics
    if os.path.exists("summary_eda_outputs/phase1_summary_table.csv"):
        phase1_df = load_and_standardize("summary_eda_outputs/phase1_summary_table.csv")
        phase1_df.drop(columns=["num_substitutions"], errors="ignore", inplace=True)
        base_df = pd.merge(base_df, phase1_df, on="peptide", how="left", suffixes=("", "_p1"))

    # Merge Step 2 specific cytokine/tier annotations
    if os.path.exists("selection_outputs/final_af3_2aa_candidates.csv"):
        af3_2aa_df = load_and_standardize("selection_outputs/final_af3_2aa_candidates.csv")
        cols_to_use = [c for c in af3_2aa_df.columns if c not in base_df.columns or c == "peptide"]
        base_df = pd.merge(base_df, af3_2aa_df[cols_to_use], on="peptide", how="left")

    # 4. Merge NetMHCpan & NetMHCstabpan Predictions (Corrected paths)
    netmhc1 = load_and_standardize("predictions/results_step1/netmhcpan/netmhcpan_predictions_full_table.csv")
    netmhc2 = load_and_standardize("predictions/results_step2/netmhcpan/netmhcpan_predictions_full_table.csv")
    netmhc_all = pd.concat([netmhc1, netmhc2], ignore_index=True).drop_duplicates(subset=["peptide"])

    stab1 = load_and_standardize("predictions/results_step1/netmhcpan/netmhcpanstab_predictions_full_table.csv")
    stab2 = load_and_standardize("predictions/results_step2/netmhcpan/netmhcpanstab_predictions_full_table.csv")
    stab_all = pd.concat([stab1, stab2], ignore_index=True).drop_duplicates(subset=["peptide"])

    neo1 = load_and_standardize("predictions/results_step1/neo_score/neo_scores.csv")
    neo2 = load_and_standardize("predictions/results_step2/neo_score/neo_scores.csv")
    neo_all = pd.concat([neo1, neo2], ignore_index=True).drop_duplicates(subset=["peptide"])

    cols_netmhc = [
        "peptide",
        "HLA-A*02:01_%Rank_EL",
        "HLA-A*02:01_Aff(nM)",
        "HLA-B*07:02_%Rank_EL",
        "HLA-C*04:01_%Rank_EL",
    ]
    cols_netmhc = [c for c in cols_netmhc if c in netmhc_all.columns]
    if cols_netmhc:
        base_df = pd.merge(base_df, netmhc_all[cols_netmhc], on="peptide", how="left")

    cols_stab = [
        "peptide",
        "HLA-A*02:01_Thalf(h)",
        "HLA-A*02:01_%Rank_Stab",
    ]
    cols_stab = [c for c in cols_stab if c in stab_all.columns]
    if cols_stab:
        base_df = pd.merge(base_df, stab_all[cols_stab], on="peptide", how="left")

    cols_neo = ["peptide", "Foreignness", "Dissimilarity"]
    cols_neo = [c for c in cols_neo if c in neo_all.columns]
    if cols_neo:
        base_df = pd.merge(base_df, neo_all[cols_neo], on="peptide", how="left")

    # 5. Merge AF3 Structural Scores
    if os.path.exists("af3/selection_outputs/af3_compiled_results.csv"):
        af3_res = load_and_standardize("af3/selection_outputs/af3_compiled_results.csv")
        cols_af3 = [c for c in ["peptide", "iptm", "ptm", "ranking_score"] if c in af3_res.columns]
        base_df = pd.merge(
            base_df,
            af3_res[cols_af3].drop_duplicates(subset=["peptide"]),
            on="peptide",
            how="left",
        )

    # 6. Merge Rosetta Relax & Rosetta Interface Analyzer
    if os.path.exists("rosetta/rosetta_candidate_summary.csv"):
        ros_summary = load_and_standardize("rosetta/rosetta_candidate_summary.csv")
        cols_ros = [
            c
            for c in [
                "peptide",
                "best_total_score",
                "mean_total_score",
                "best_fa_rep",
                "best_hbond_sc",
            ]
            if c in ros_summary.columns
        ]
        base_df = pd.merge(
            base_df,
            ros_summary[cols_ros].rename(
                columns={
                    "best_total_score": "rosetta_relax_best_total",
                    "mean_total_score": "rosetta_relax_mean_total",
                }
            ),
            on="peptide",
            how="left",
        )

    ros_int_summary = parse_rosetta_interface()
    if not ros_int_summary.empty:
        base_df = pd.merge(base_df, ros_int_summary, on="peptide", how="left")

    # Deduplicate columns if any duplicate headers exist
    base_df = base_df.loc[:, ~base_df.columns.duplicated()]

    # Save Master Table
    out_file = "master_34_candidates_annotated.csv"
    base_df.to_csv(out_file, index=False)
    print(
        f"Successfully created {out_file} with shape {base_df.shape} across {len(base_df.columns)} columns."
    )
    return base_df


if __name__ == "__main__":
    df_master = build_master_table()
