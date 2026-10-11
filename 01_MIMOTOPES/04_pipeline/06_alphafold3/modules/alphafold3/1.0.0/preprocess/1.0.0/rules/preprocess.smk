import scripts.utils as utils
import pandas as pd
from scripts.utils.path_builder import *
writer = Writer(config["benchmark"], directory = "results", pipeline = "alphafold3", module = "preprocess")

rulename = "process_peptides"
rule process_peptides:
    """Process Peptide."""
    version: "1.0.0"
    # message: "Convert the peptide {wildcards.sample} into a Fasta File."
    # shadow: "full"
    input:
        json = lambda wildcards: str(config["samples"].loc[wildcards.sample, "file"])
    output:
        fasta = writer.outputs("{sample}", "{threads}", rulename, "{sample}/{sample}_data.json")
    params:
        script = "/rsrch6/home/hema_bio-Malignan/adp_shukla/pipeline_resources/alphafold3/alphafold3/run_alphafold.py"
    # group: "alphafold3_group"
    threads: lambda wildcards, rulename=rulename: int(wildcards.threads) if config["benchmark"] else config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = writer.io(rulename).log_stdout,
        stderr = writer.io(rulename).log_stderr
    benchmark:
        writer.io(rulename).benchmark
    # singularity: config["containers"]["alphafold3"]
    shell:
        """
        apptainer exec --net --network none --nv \
        -B /rsrch4:/rsrch4 \
        -B /rsrch6:/rsrch6 \
        -B /rsrch9:/rsrch9 \
        -B "/rsrch6/home/hema_bio-Malignan/adp_shukla/pipeline_resources/alphafold3/data/params:/data/models" \
        -B "/rsrch6/home/hema_bio-Malignan/adp_shukla/pipeline_resources/alphafold3/public_databases:/data/public_databases" \
        /rsrch6/home/hema_bio-Malignan/adp_shukla/pipeline_resources/alphafold3/alphafold3.sif \
        python {params.script} \
        --json_path={input.json} \
        --model_dir=/data/models \
        --db_dir=/data/public_databases \
        --output_dir=results/{wildcards.sample}/alphafold3/preprocess/process_peptides/
        """
