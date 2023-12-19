"""
All target download files are declared here
"""


def list_assembly_downloads(wildcards):
    return list(
        pd.read_csv(
            checkpoints.download_assemblies.get(**wildcards).output[0], sep="\t"
        )["path"]
    )


TargetDownloads = [
    os.path.join(dir.out.base, "download/assembly_summary.tsv"),
    list_assembly_downloads,
]
