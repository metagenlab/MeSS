"""
All target download files are declared here
"""

TargetDownloads = [
    os.path.join(dir.out.base, "uniq_entries.tsv"),
    os.path.join(dir.out.base, "assembly_finder/assembly_summary.tsv"),
    os.path.join(dir.out.base, "assembly_finder/taxonomy.tsv"),
]
