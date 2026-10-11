import time
import pandas as pd
from pathlib import Path
import sys

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PIPELINE_ROOT))

from lib.logging_utils import setup_logger
from lib.restriction_parser import parse_restrictions
from lib.viz_utils import generate_restriction_heatmap
from generate_variants import MimotopeGenerator

def main():
    log_file = PIPELINE_ROOT / "logs" / "02_mimotopes_generation.log"
    logger = setup_logger("02_mimotopes_generation", log_file=log_file)
    logger.info("Starting Module 02: In Silico Mimotope Generation")

    input_dir = PIPELINE_ROOT / "01_input"
    output_dir = PIPELINE_ROOT / "02_mimotopes_generation" / "output"
    figures_dir = PIPELINE_ROOT / "02_mimotopes_generation" / "figures"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    peptides_csv = input_dir / "peptides.csv"
    peptides_df = pd.read_csv(peptides_csv)

    for idx, row in peptides_df.iterrows():
        wt_seq = row["peptide"].strip()
        restr_filename = str(row.get("restrictions", "")).strip()
        
        logger.info(f"Processing Peptide [{idx+1}/{len(peptides_df)}]: {wt_seq}")
        start_time = time.time()

        constraints = None
        if restr_filename and restr_filename != "nan":
            restr_path = input_dir / restr_filename
            if restr_path.exists():
                constraints = parse_restrictions(restr_path, wt_seq)
                logger.info(f"Parsed constraints for all {len(constraints)} positions from {restr_filename}")

                heatmap_path = figures_dir / f"{wt_seq}_restriction_heatmap.png"
                generate_restriction_heatmap(wt_seq, constraints, heatmap_path)
                logger.info(f"Saved corrected restriction heatmap: {heatmap_path}")

        # Generation
        generator = MimotopeGenerator(wt_seq)
        variants_df = generator.generate_all(constraints=constraints)

        # Verification Check
        if constraints:
            generator.verify_ps_cpl_constraints(variants_df, constraints, wt_seq)
            logger.info("VERIFICATION PASSED: 100% of PS-CPL variants strictly respect positional constraints.")

        # Save Output
        out_csv = output_dir / f"{wt_seq}_in_silico_variants.csv"
        variants_df.to_csv(out_csv, index=False)
        
        elapsed = time.time() - start_time
        logger.info(f"Completed {wt_seq} in {elapsed:.2f}s. Saved {len(variants_df):,} variants to {out_csv}")

    logger.info("Module 02 execution completed successfully.")

if __name__ == "__main__":
    main()
