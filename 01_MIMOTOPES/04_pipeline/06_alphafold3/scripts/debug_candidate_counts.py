# scripts/debug_candidate_counts.py
from pathlib import Path
import pandas as pd

def debug_counts():
    meta_path = Path("../05_in_silico_variants_selection/output/in_silico_variants_selected.csv")
    af3_csv = Path("output/af3_compiled_results.csv")
    plddt_csv = Path("output/peptide_anchor_stability.csv")
    results_dir = Path("results")

    print("=== 1. Metadata Target Counts ===")
    if meta_path.exists():
        meta_df = pd.read_csv(meta_path)
        print(meta_df["category"].value_counts())
    else:
        print(f"Missing: {meta_path}")

    print("\n=== 2. AF3 Output Folders on Disk ===")
    cif_files = list(results_dir.glob("**/*_model.cif"))
    print(f"Total CIF models found: {len(cif_files)}")

    print("\n=== 3. Compiled AF3 Results per Category ===")
    if af3_csv.exists():
        comp_df = pd.read_csv(af3_csv)
        print(comp_df["category"].value_counts())
    else:
        print(f"Missing: {af3_csv}")

    print("\n=== 4. Unique Peptides in pLDDT Stability Output ===")
    if plddt_csv.exists():
        plddt_df = pd.read_csv(plddt_csv)
        unique_peps = plddt_df.groupby("category")["peptide"].nunique()
        print(unique_peps)
    else:
        print(f"Missing: {plddt_csv}")

if __name__ == "__main__":
    debug_counts()
