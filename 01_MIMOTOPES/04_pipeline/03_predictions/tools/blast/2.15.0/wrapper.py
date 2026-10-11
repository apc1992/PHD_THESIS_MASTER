from snakemake.shell import shell
import os
from datetime import datetime

# Get variables from snakemake
input = snakemake.input[0]
output = snakemake.output[0]
threads = snakemake.threads
log = snakemake.log[0]
container = snakemake.params.get("container")
task = snakemake.params.get("task")
database = snakemake.params.get("database")
output_format = snakemake.params.get("output_format")
args = snakemake.params.get("args")

# Ensure the output directory exists
os.makedirs(os.path.dirname(output), exist_ok=True)

# Shell command
blast_cmd = f"""
singularity exec --net --network none -B $(pwd) -B /rsrch4:/rsrch4 -B /rsrch5:/rsrch5 -B /rsrch6:/rsrch6 -B /rsrch9:/rsrch9 {container} \
blastp {task} -db {database} -query {input} -outfmt {output_format} {args} -num_threads {threads} > {output}
"""

shell(blast_cmd)