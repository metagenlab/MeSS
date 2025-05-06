# config file
import attrmap as ap


configfile: os.path.join(workflow.basedir, "config", "config.yaml")


config = ap.AttrMap(config)


# directories
include: "directories.smk"
# functions
include: "functions.smk"


onsuccess:
    copy_log_file()


onerror:
    copy_log_file()


# common options
INPUT = os.path.abspath(str(config.args.input))
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "mess.log")
THREADS = config.args.threads
TAXONKIT = config.args.taxonkit


include: "setup.smk"
