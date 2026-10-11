import os
from snakemake.shell import shell

# Get variables from snakemake
input_peps = snakemake.input.get("peps")
human_alignment = snakemake.input.get("human_alignment")
iedb_alignment = snakemake.input.get("iedb_alignment")
output = snakemake.output.get("pep_scores")
container = snakemake.params.get("container")

# Get the absolute path to neo_scores.py sitting next to this wrapper
script_dir = os.path.dirname(__file__)
script_path = os.path.join(script_dir, "neo_scores.py")

# Construct the command 
shell_command = f"""  
    singularity exec --net --network none -B $(pwd) -B /rsrch4:/rsrch4 -B /rsrch5:/rsrch5 -B /rsrch6:/rsrch6 -B /rsrch9:/rsrch9 {container} \
    python {script_path} --input {input_peps} --human_hits {human_alignment} --iedb_hits {iedb_alignment} --output {output}
    """

# Execute the combined command
shell(shell_command)