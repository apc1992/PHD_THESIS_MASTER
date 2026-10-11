import os

## Module config files: Rule related references and references
configfile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "template", "1.0.0", "config", "config.yaml")
configfile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "netmhcpan", "4.2.0", "config", "config.yaml")
configfile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "netmhcstab", "1.0.1", "config", "config.yaml")
configfile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "antigen_scoring", "1.0.0", "config", "config.yaml")

## Load modules

module neoscoreT:
    snakefile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "template", "1.0.0", "Snakefile")
    config: config

module netmhcpan:
    snakefile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "netmhcpan", "4.2.0", "Snakefile")
    config: config

module netmhcstab:
    snakefile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "netmhcstab", "1.0.1", "Snakefile")
    config: config

module antigen_scoring:
    snakefile: os.path.join(os.getcwd(), "modules", "neoscore", "1.0.0", "antigen_scoring", "1.0.0", "Snakefile")
    config: config

# Rules from the modules
use rule * from netmhcpan as netmhcpan_*
use rule * from netmhcstab as netmhcstab_*
use rule * from antigen_scoring as antigen_scoring_*
