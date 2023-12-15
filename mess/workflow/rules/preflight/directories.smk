"""
Ensure consistent directory names

A lot inspired from hybracter and hecatomb
"""

import attrmap as ap

dir = ap.AttrMap()

# Workflow dirs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

# Output location
try:
    assert (ap.utils.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "mess_out"

# logs and benchmarks
dir.out.bench = os.path.join(dir.out.base, "benchmarks")
dir.out.stderr = os.path.join(dir.out.base, "stderr")
