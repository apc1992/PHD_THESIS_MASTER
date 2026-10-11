import itertools
import pandas as pd
from typing import Dict, List, Set, Tuple

STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

class MimotopeGenerator:
    def __init__(self, wt_seq: str):
        self.wt_seq = wt_seq
        self.seq_len = len(wt_seq)

    def generate_all(self, constraints: Dict[int, List[str]] = None) -> pd.DataFrame:
        records: List[Dict] = []
        seen_seqs: Set[str] = set()

        # 0. Wild-Type Baseline
        records.append(self._create_record(self.wt_seq, [], "WT"))
        seen_seqs.add(self.wt_seq)

        # 1. 1-AA Substitutions
        for pos in range(self.seq_len):
            wt_aa = self.wt_seq[pos]
            for aa in STANDARD_AAS:
                if aa == wt_aa:
                    continue
                mimo = list(self.wt_seq)
                mimo[pos] = aa
                mimo_str = "".join(mimo)
                
                if mimo_str not in seen_seqs:
                    seen_seqs.add(mimo_str)
                    records.append(self._create_record(mimo_str, [(pos, wt_aa, aa)], "1AA"))

        # 2. 2-AA Substitutions
        for pos1, pos2 in itertools.combinations(range(self.seq_len), 2):
            wt_aa1, wt_aa2 = self.wt_seq[pos1], self.wt_seq[pos2]
            for aa1 in STANDARD_AAS:
                if aa1 == wt_aa1:
                    continue
                for aa2 in STANDARD_AAS:
                    if aa2 == wt_aa2:
                        continue
                    mimo = list(self.wt_seq)
                    mimo[pos1] = aa1
                    mimo[pos2] = aa2
                    mimo_str = "".join(mimo)

                    if mimo_str not in seen_seqs:
                        seen_seqs.add(mimo_str)
                        mutations = [(pos1, wt_aa1, aa1), (pos2, wt_aa2, aa2)]
                        records.append(self._create_record(mimo_str, mutations, "2AA"))

        # 3. PS-CPL Guided Substitutions
        if constraints:
            # Filter positions that have >1 allowed residue
            variable_positions = [p for p in range(self.seq_len) if len(constraints[p]) > 1]
            aa_choices = [constraints[p] for p in variable_positions]

            for residues in itertools.product(*aa_choices):
                mimo = list(self.wt_seq)
                mutations: List[Tuple[int, str, str]] = []

                for pos, mut_aa in zip(variable_positions, residues):
                    wt_aa = self.wt_seq[pos]
                    if mut_aa != wt_aa:
                        mimo[pos] = mut_aa
                        mutations.append((pos, wt_aa, mut_aa))

                mimo_str = "".join(mimo)
                num_muts = len(mutations)

                if num_muts >= 3 and mimo_str not in seen_seqs:
                    seen_seqs.add(mimo_str)
                    records.append(self._create_record(mimo_str, mutations, "PS-CPL"))

        return pd.DataFrame(records)

    @staticmethod
    def verify_ps_cpl_constraints(df: pd.DataFrame, constraints: Dict[int, List[str]], wt_seq: str) -> bool:
        """
        Rigorously asserts that every PS-CPL generated mimotope contains ONLY 
        allowed amino acids at every position according to restrictions.
        """
        ps_cpl_df = df[df["category"] == "PS-CPL"]
        
        for idx, row in ps_cpl_df.iterrows():
            mimo = row["mimotope"]
            for pos_idx, aa in enumerate(mimo):
                allowed_set = set(constraints[pos_idx])
                if aa not in allowed_set:
                    raise AssertionError(
                        f"Constraint Violation at row {idx}! Mimotope '{mimo}' has '{aa}' "
                        f"at position {pos_idx+1}, which is not in allowed set: {allowed_set}"
                    )
        return True

    def _create_record(self, mimotope: str, mutations: List[Tuple[int, str, str]], category: str) -> Dict:
        if not mutations:
            return {
                "peptide": self.wt_seq,
                "mimotope": mimotope,
                "substitution": "WT",
                "position": "0",
                "mut_aminoacid": "WT",
                "num_substitutions": 0,
                "category": category
            }

        substitutions = [f"{wt_aa}->{mut_aa}" for _, wt_aa, mut_aa in mutations]
        positions = [str(pos + 1) for pos, _, _ in mutations]
        mut_aas = [mut_aa for _, _, mut_aa in mutations]

        return {
            "peptide": self.wt_seq,
            "mimotope": mimotope,
            "substitution": ";".join(substitutions),
            "position": ";".join(positions),
            "mut_aminoacid": ";".join(mut_aas),
            "num_substitutions": len(mutations),
            "category": category
        }
