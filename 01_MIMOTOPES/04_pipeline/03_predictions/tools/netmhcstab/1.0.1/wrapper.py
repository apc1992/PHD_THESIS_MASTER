from snakemake.shell import shell
import os

# Get variables from snakemake
input = snakemake.input[0]
output = snakemake.output[0]
threads = snakemake.threads
log = snakemake.log[0]
container = snakemake.params.get("container")
alleles = snakemake.params.get("alleles")

# Ensure the output directory exists
os.makedirs(os.path.dirname(output), exist_ok=True)

# Shell command
shell_command = f"""  
    singularity exec --net --network none -B $(pwd) -B /rsrch4:/rsrch4 -B /rsrch5:/rsrch5 -B /rsrch6:/rsrch6 -B /rsrch9:/rsrch9 {container} \
    /usr/local/bin/netMHCstabpan-1.0/netMHCstabpan -f "{input}" -inptype 1 -a "{alleles}" > "{output}"
"""

# Execute the combined command
shell(shell_command)