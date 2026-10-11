from pathlib import Path

# Target Wild-Type Sequence & Experimental Metrics
WT_SEQ = "QCFFCFKEL"
CYTOKINES = ['CD4_41BB', 'CD4_IFNg', 'CD8_41BB', 'CD8_IFNg']
DONORS = ['donor_1', 'donor_2', 'donor_3']
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
POSITIONS = list(range(1, len(WT_SEQ) + 1))
EPSILON = 1e-5

# Global Pipeline Paths
ROOT_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT_DIR / "01_input"
PRED_DIR = ROOT_DIR / "03_predictions/results"

# Visualization Styling Defaults
TITLE_SIZE = 18
AXIS_LABEL_SIZE = 14
LEGEND_SIZE = 12
TICK_LABEL_SIZE = 11
