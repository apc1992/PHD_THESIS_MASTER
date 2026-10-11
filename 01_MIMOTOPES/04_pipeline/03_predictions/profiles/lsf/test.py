from CookieCutter import CookieCutter

# Get the raw string from the profile's config.yaml:
cluster_cfg_path = CookieCutter.cluster_config_path()
print()
print(cluster_cfg_path)
print()
# e.g., "pipelines/IsoSplit/1.0.0/profiles/lsf/cluster.yaml"

# Or load the referenced cluster YAML file into a dict:
cluster_cfg = CookieCutter.cluster_config_yaml()
print()
print(cluster_cfg)
print()
# e.g., {"__default__": {...}, "rule_a": {...}, ...}
