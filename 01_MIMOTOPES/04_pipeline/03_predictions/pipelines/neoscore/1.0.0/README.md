# MHC-I predictions
![NetMHCPan](https://img.shields.io/badge/NetMHCPan-v4.1b-blue)
![MHC-I Predictions](https://img.shields.io/badge/MHC_I_Predictions-v.1.0.0-8A2BE2)

Snakemake workflow to run MHC-I predictions for a set of peptides, and summarize the Binding Affinity (BA) / Eluted Ligand (EL) scores and ranks.

## Workflow

<p align="center"><img src="resources/workflow.png" alt="Workflow" width="800"></p>

## Usage

Ensure you have `Snakemake` and `Singularity` installed, as this pipeline uses containerized software. Modify the `config.yaml` to reflect your system configuration and desired parameters.

Bear in mind that you will need to have an `LSF profile` and `Snakemake` installed in your environment to run the workflow.

## 1. Creating the LSF profile

```
# create configuration directory that snakemake searches for profiles
profile_dir="${HOME}/.config/snakemake"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/lsf"
cookiecutter --output-dir "$profile_dir" "$template"
```
Different questions will appear, we just highlight the important ones here:

```
LSF_UNIT_FOR_LIMITS: 3 (GB)
jobs: 10000
default_queue: e80medium
max_jobs_per_second [10]: 100
other options by default
```

## 2. Run on Seadragon
1. Clone the repository in your working environment

`git clone git@github.mdanderson.org:shukla-lab/MHC-I-Predictions.git`

`cd MHC-I-Predictions`

2. Activate your environment with snakemake, e.g. through conda

`conda activate snakemake`

or append the commands in your submission LSF script

`source ~/.bashrc`

`conda activate snakemake`

3. Add the path to your input file in the configuration file `config["input"]`. Optionally, check other default parameters.

4. Submit the LSF script to the job scheduler

`bsub < run.lsf` 

## Input libraries as a CSV file

The workflow is designed to receive the output from the [Fragpipe workflow](https://github.mdanderson.org/shukla-lab/fragpipe/tree/main), but it accepts any tab-separated table one column named either `peptide` or `pep`, case-insensitive, that contains the peptides to be queried. The peptides shorter than 8 or longer than 12 amino acids will be filtered out. 

## What's pending?

None

## Contributions

We welcome contributions to this project! If you would like to contribute, please follow these steps:

Fork the repository.
Create a new branch for your feature or bug fix.
Make your changes and commit them with descriptive messages.
Push your changes to your forked repository.
Submit a pull request, explaining the changes you have made.
Please ensure that your contributions adhere to our code of conduct.

## Contact 

CYu5@mdanderson.org

JJGallegos@mdanderson.org

BPan2@mdanderson.org

APeris@mdanderson.org