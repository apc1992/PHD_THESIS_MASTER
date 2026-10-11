import sys
from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib.logging_utils import setup_logger
from cytokine_scorer import validate_and_extract_cytokine_cols, compute_wt_relative_ratios, score_mimotopes
from heatmap_builder import generate_donor_heatmaps
from publication_plots import generate_binding_summary_figure, generate_functional_immunogenicity_figure

def run_data_summary():
    base_dir = Path(__file__).resolve().parents[1]
    root_dir = base_dir.parent
    logger = setup_logger("04_data_summary", root_dir / "logs" / "04_data_summary.log")
    
    logger.info("Initializing step 04_data_summary pipeline module...")

    peptides_config_path = root_dir / "01_input" / "peptides.csv"
    if not peptides_config_path.exists():
        raise FileNotFoundError(f"Missing master input configuration: {peptides_config_path}")

    peptides_df = pd.read_csv(peptides_config_path)

    for _, p_row in peptides_df.iterrows():
        peptide_seq = str(p_row['peptide']).strip()
        cytokines_file = str(p_row['cytokines']).strip()
        
        logger.info(f"Processing candidate peptide: {peptide_seq}")

        cytokines_path = root_dir / "01_input" / cytokines_file
        variants_path = root_dir / "02_mimotopes_generation" / "output" / f"{peptide_seq}_in_silico_variants.csv"
        neo_path = root_dir / "03_predictions" / "results" / "neo_score" / "neo_scores.csv"
        netmhc_path = root_dir / "03_predictions" / "results" / "netmhcpan" / "netmhcpan_predictions_full_table.csv"
        stab_path = root_dir / "03_predictions" / "results" / "netmhcpan" / "netmhcpanstab_predictions_full_table.csv"

        if not variants_path.exists():
            logger.error(f"Expected variant output missing: {variants_path}")
            continue

        # 1. Process Cytokine Activity Ratios & Heatmaps
        if cytokines_path.exists():
            cytokines_df = pd.read_csv(cytokines_path)
            cyto_cols = validate_and_extract_cytokine_cols(cytokines_df)
            
            logger.info(f"Building relative activity heatmaps for {peptide_seq}...")
            generate_donor_heatmaps(cytokines_df, peptide_seq, cyto_cols, base_dir / "figures", peptide_seq)
            
            logger.info(f"Calculating combinatorial activity ratios for {peptide_seq}...")
            ratios = compute_wt_relative_ratios(cytokines_df, peptide_seq, cyto_cols)
        else:
            logger.warning(f"Cytokines input missing at {cytokines_path}. Skipping cytokine scoring.")
            ratios, cyto_cols = {}, []

        # 2. Annotate Variants with Predictions and Functional Scores
        logger.info(f"Loading variants and merging predictions for {peptide_seq}...")
        variants_df = pd.read_csv(variants_path)

        if ratios:
            scored_df = score_mimotopes(variants_df, peptide_seq, ratios, cyto_cols)
            variants_df = pd.merge(variants_df, scored_df, on="mimotope", how="left")

        # Rename 'Peptide' to 'mimotope' in prediction tables before merging to prevent duplicate Peptide_x/Peptide_y columns
        for pred_path in [neo_path, netmhc_path, stab_path]:
            if pred_path.exists():
                pred_df = pd.read_csv(pred_path)
                if 'Peptide' in pred_df.columns:
                    pred_df = pred_df.rename(columns={'Peptide': 'mimotope'})
                
                cols_to_use = [c for c in pred_df.columns if c not in variants_df.columns or c == 'mimotope']
                variants_df = pd.merge(variants_df, pred_df[cols_to_use], on='mimotope', how='left')

        out_csv = base_dir / "output" / f"{peptide_seq}_in_silico_variants_metrics.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        variants_df.to_csv(out_csv, index=False)
        logger.info(f"Exported master metrics file: {out_csv}")

        # 3. Generate Descriptive Figures
        logger.info(f"Generating summary figures for {peptide_seq}...")
        generate_binding_summary_figure(variants_df, base_dir / "figures", peptide_seq)
        generate_functional_immunogenicity_figure(variants_df, base_dir / "figures", peptide_seq)

    logger.info("Module 04_data_summary execution completed successfully.")

if __name__ == "__main__":
    run_data_summary()
