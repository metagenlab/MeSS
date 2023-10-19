import gzip
from Bio import SeqIO
import pandas as pd
from collections import Counter
seqtech = snakemake.params[0]


def get_read_counts(fastq_gz, ss):
    """
    function that returns a table with each genome read count
    """
    with gzip.open(fastq_gz, 'rt') as handle:
        fastq_parser = SeqIO.parse(handle, 'fastq')
        headers = [record.description for record in fastq_parser]
        if ss == 'illumina':
            headers = [header.split('-')[-1].split('.')[0] for header in headers]
        else:
            headers = [header.split('.')[0] for header in headers]
        c = dict(Counter(headers))
        rc = pd.DataFrame.from_dict(c, orient='index', columns=['true_read_counts'])
        rc.index.set_names('Assembly', inplace=True)
        rc.reset_index(inplace=True)
    return rc


read_counts = get_read_counts(snakemake.input.reads, seqtech)
summary = pd.read_csv(snakemake.input.table, sep='\t')
summary['Assembly'] = [assemblyname.split('.')[0] for assemblyname in summary['AssemblyNames']]
merged = pd.merge(summary, read_counts, on='Assembly')
cols = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'Taxid', 'true_read_counts']
merged[cols].to_csv(snakemake.output[0], sep='\t', index=None)
krona_cols = ['true_read_counts', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
krona = merged[krona_cols]
krona.set_index('true_read_counts').to_csv(snakemake.output.krona, sep='\t')
