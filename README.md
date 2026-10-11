# PHD_MASTER — AI WORKSPACE SETUP

This is a **non-destructive starter package** for the AI-assisted PhD workflow.

It does not contain or modify your existing research data, manuscripts, code, figures or thesis files. Copy the `00_AI` contents and run `SETUP_EXISTING_PHD_MASTER.sh` from the root of your existing `PHD_MASTER` directory.

## Final scientific priority

1. Mimotopes: full computational/project review and remediation; prepare Thesis Part III while final immunogenicity experiments are pending.
2. Pan-cancer: complete the minimum validated analysis needed for Thesis Parts I–II; build initial thesis figures/tables. The full Publication 2 manuscript is deferred.
3. Thesis: write and review the complete book-style dissertation until submission-ready.
4. After thesis submission: return to publication finalization.

## Thesis structure

1. Current approaches for antigen discovery in cancer
2. Application of a cutting-edge antigen discovery platform in a large pan-cancer cohort
3. Optimization of antigen candidates with mimotopes

## How to deploy

From your existing `PHD_MASTER` directory:

```bash
# 1. Back up the current directory first if it contains important unpublished work.

# 2. Copy the prepared 00_AI directory into PHD_MASTER, reviewing conflicts.
#    Do not blindly overwrite existing files you have modified since the last setup.

# 3. Copy the setup script into PHD_MASTER and run:
chmod +x SETUP_EXISTING_PHD_MASTER.sh
./SETUP_EXISTING_PHD_MASTER.sh
```

The script is deliberately non-destructive. It creates missing directories and copies the currently known PanCan root-level files into their intended folders only when a target file does not already exist.

## First working sequence

### Pan-cancer inventory

Open:

`03_PANCAN/cohort/PANCAN_SAMPLE_ANALYSIS_MATRIX.xlsx`

Populate one sample per row, replace AN01–AN30 with the real analyses, and define each analysis in `Analysis_Register`.

### Mimotopes

Start with a **read-only audit** of:

`01_MIMOTOPES/04_pipeline`

Then map code/output → manuscript → figures/tables → thesis Part III.

### Thesis

Start with the institutional rules in:

`02_PHD_THESIS/rules`

Then build the thesis outline and evidence map before doing a full stylistic rewrite.

## Local-vs-HPC workflow

Use the Mac for writing, literature work, lightweight coding, figure generation and project organization. Keep heavy computation on HPC and export only the appropriate derived outputs needed for downstream writing/figure work.

Do not connect the AI workflow to institutional HPC until the relevant institutional permissions and technical constraints have been clarified.
