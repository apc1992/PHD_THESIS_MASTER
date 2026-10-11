import os

# Load module configurations
configfile: os.path.join(os.getcwd(), "pipelines", "alphafold3", "1.0.0", "config", "config.yaml")
configfile: os.path.join(os.getcwd(), "modules", "alphafold3", "1.0.0", "preprocess", "1.0.0", "config", "config.yaml")

# Load modules 
module preprocess:
    snakefile: os.path.join(os.getcwd(), "modules", "alphafold3", "1.0.0", "preprocess", "1.0.0", "Snakefile")
    config: config

# Rules from the modules
use rule * from preprocess as preprocess_*