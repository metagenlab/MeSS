import os
include: 'rules/find_assemblies.rules'
community_name=config['community_name']
def downloaded_list(wildcards):
    checkpoint_output = checkpoints.download_assemblies.get(**wildcards).output[0]
    directory = '/'.join((checkpoint_output.split('/')[0:2]))
    return expand(f'assembly_gz/{community_name}/{{i}}_genomic.fna.gz',
           i=glob_wildcards(os.path.join(directory, '{i}_genomic.fna.gz')).i)

rule all_download:
    input: f"{community_name}-assemblies-summary.tsv",
            downloaded_list,
            f"assembly_gz/{community_name}/{community_name}.done"