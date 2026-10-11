import pandas as pd
import numpy as np
import argparse
import os
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load

def get_sw_score(peptide, hit_seq):
    """
    Calculates Smith-Waterman local alignment score between a peptide 
    and a database hit.
    """
    # Clean sequences: Biopython .score() needs raw AAs, not alignment strings with gaps
    p = str(peptide).upper().replace('-', '').replace('*', '').strip()
    h = str(hit_seq).upper().replace('-', '').replace('*', '').strip()
    
    if not p or not h:
        return 0

    # Alphabet check: BLOSUM62 does not contain U, X, B, Z, etc.
    # Non-standard AAs will trigger a ValueError in Biopython.
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    if not set(p).issubset(valid_aa) or not set(h).issubset(valid_aa):
        return 0

    # Setup Aligner with Łuksza et al. 2017 parameters
    aligner = PairwiseAligner()
    aligner.substitution_matrix = load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = 'local' 
    
    try:
        return aligner.score(p, h)
    except ValueError:
        return 0

def calculate_partition(scores, a, k=4.86936):
    """
    Implements the fitness model partition function: 
    Z = sum(exp[-k(a - s)]) / (1 + sum(exp[-k(a - s)]))
    """
    if not scores: 
        return 0.0
    
    scores = np.array(scores)
    # The exponent: -k * (threshold - score)
    be = -k * (a - scores)
    sum_exp = np.sum(np.exp(be))
    return sum_exp / (1 + sum_exp)

def load_hits(hit_file):
    """
    Loads BLAST output and cleans subject sequences.
    """
    hits = {}
    if not os.path.exists(hit_file) or os.path.getsize(hit_file) == 0:
        return hits
        
    try:
        # Expected format: qseqid, sseqid, qseq, sseq
        df = pd.read_csv(hit_file, names=['q_id', 's_id', 'q_seq', 's_seq'], sep=',')
        for _, row in df.iterrows():
            try:
                idx = int(row['q_id'])
                if idx not in hits: 
                    hits[idx] = []
                
                # We store the cleaned s_seq. 
                # We don't need the alignment gaps because SW will re-align.
                s_clean = str(row['s_seq']).upper().replace('-', '').replace('*', '').strip()
                if s_clean:
                    hits[idx].append(s_clean)
            except (ValueError, TypeError):
                continue
    except Exception as e:
        print(f"Warning: Issue parsing {hit_file}: {e}")
        
    return hits

def main():
    parser = argparse.ArgumentParser(description="Calculate Neoantigen Fitness Scores")
    parser.add_argument("--input", required=True, help="Input peptide table (txt/csv)")
    parser.add_argument("--human_hits", required=True, help="BLAST hits against human proteome")
    parser.add_argument("--iedb_hits", required=True, help="BLAST hits against IEDB")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    # Load input peptides
    df = pd.read_csv(args.input, names=['Peptide'])
    peptides = df['Peptide'].tolist()
    
    # Load alignments
    h_hits = load_hits(args.human_hits)
    i_hits = load_hits(args.iedb_hits)

    # 1. Foreignness calculation (IEDB hits, threshold a=26)
    # Rationale: Measures mimicry of known pathogens
    foreignness = []
    for i, p in enumerate(peptides):
        scores = [get_sw_score(p, h) for h in i_hits.get(i, [])]
        foreignness.append(calculate_partition(scores, 26))
    df['Foreignness'] = foreignness
    
    # 2. Dissimilarity calculation (Human hits, threshold a=32)
    # Rationale: Measures how different the peptide is from the self-proteome
    dissimilarity = []
    for i, p in enumerate(peptides):
        scores = [get_sw_score(p, h) for h in h_hits.get(i, [])]
        # 1 - Z(human) because higher similarity to self means LOWER fitness
        dissimilarity.append(1 - calculate_partition(scores, 32))
    df['Dissimilarity'] = dissimilarity

    # Final Output
    df.to_csv(args.output, index=False)
    print(f"Successfully processed {len(df)} peptides.")
    print(f"Results saved to: {args.output}")

if __name__ == "__main__":
    main()