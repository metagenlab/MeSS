"""
Ensure consistent directory names

A lot inspired from hybracter and hecatomb
"""

import attrmap as ap
import yaml

dir = ap.AttrMap()

# Env paths
dir.conda = os.path.join(workflow.basedir, "envs", "conda")
containers = ap.AttrMap(
    yaml.safe_load(open(os.path.join(workflow.basedir, "envs", "containers.yml"), "r"))
)
# Workflow dirs
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

# Base output location

try:
    assert (ap.utils.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "results"

# Outdirs
dir.out.processing = os.path.join(dir.out.base, "processing")
dir.out.cat = os.path.join(dir.out.base, "cat")
dir.out.shuffle = os.path.join(dir.out.base, "shuffle")
dir.out.short = os.path.join(dir.out.processing, "short")
dir.out.long = os.path.join(dir.out.processing, "long")
dir.out.fastq = os.path.join(dir.out.base, "fastq")
dir.out.bam = os.path.join(dir.out.base, "bam")
dir.out.tax = os.path.join(dir.out.base, "tax")


# Logs versions and benchmarks
dir.out.bench = os.path.join(dir.out.base, "benchmarks")
dir.out.logs = os.path.join(dir.out.base, "logs")
dir.out.versions = os.path.join(dir.out.base, "versions")
