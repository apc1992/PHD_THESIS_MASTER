import logging
from pathlib import Path
import pandas as pd
from Bio.PDB import MMCIFParser

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def extract_chain_plddt(cif_path: Path, chain_id: str = "C") -> list:
    """Extracts pLDDT scores per residue for a given chain from mmCIF."""
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("af3_model", str(cif_path))
    except Exception as e:
        logging.error(f"Failed parsing {cif_path}: {e}")
        return []

    model = structure[0]
    if chain_id not in model:
        return []

    residues = []
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        res_num = residue.id[1]
        res_name = residue.get_resname()
        atom = residue["CA"] if "CA" in residue else list(residue.get_atoms())[0]
        plddt = atom.get_bfactor()

        residues.append({"position": res_num, "amino_acid": res_name, "plddt": plddt})
    return residues

def main(
    results_dir: Path = Path("results"),
    compiled_csv: Path = Path("output/af3_compiled_results.csv"),
    output_csv: Path = Path("output/peptide_anchor_stability.csv"),
):
    logging.info("Extracting peptide pLDDT stability profiles...")
    cif_files = list(results_dir.glob("**/*_model.cif"))
    logging.info(f"Found {len(cif_files)} CIF model files.")

    records = []
    for cif_path in cif_files:
        folder_name = cif_path.parent.name
        parts = folder_name.split("_")
        peptide_seq = parts[-1].strip().upper()

        residues = extract_chain_plddt(cif_path, chain_id="C")
        seq_len = len(residues)

        for res in residues:
            pos = res["position"]
            anchor_type = "P2_Anchor" if pos == 2 else ("P-Omega_Anchor" if pos == seq_len else "Non_Anchor")
            records.append({
                "job_name": folder_name,
                "peptide": peptide_seq,
                "length": seq_len,
                "residue_position": f"P{pos}",
                "amino_acid": res["amino_acid"],
                "plddt": res["plddt"],
                "anchor_type": anchor_type,
            })

    plddt_df = pd.DataFrame(records)
    if plddt_df.empty:
        logging.warning("No residue pLDDT profiles extracted.")
        return

    # Attach metadata by joining on job_name first, falling back to peptide sequence
    if compiled_csv.exists():
        comp_df = pd.read_csv(compiled_csv)
        meta_cols = [c for c in ["job_name", "category", "Selection_Tier", "mimotope", "Peptide"] if c in comp_df.columns]
        
        if "job_name" in comp_df.columns and "job_name" in plddt_df.columns:
            meta_df = comp_df[meta_cols].drop_duplicates(subset=["job_name"])
            plddt_df = pd.merge(plddt_df, meta_df, on="job_name", how="left")
        else:
            match_col = "mimotope" if "mimotope" in comp_df.columns else "Peptide"
            meta_df = comp_df[meta_cols].drop_duplicates(subset=[match_col, "category"])
            plddt_df = pd.merge(plddt_df, meta_df, left_on="peptide", right_on=match_col, how="left")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    plddt_df.to_csv(output_csv, index=False)
    logging.info(f"Exported {len(plddt_df)} residue stability records to: {output_csv}")

if __name__ == "__main__":
    main()
