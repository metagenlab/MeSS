import pandas as pd
import numpy as np
import logging
import random
replicate = snakemake.wildcards.rep
seed = snakemake.params.seed
random.seed(seed)
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
inputval = snakemake.params.value
if inputval == 'Reads' or inputval == 'Coverage':
    totalreads = 0
else:
    totalreads = snakemake.config['total_reads']
sdr = snakemake.config["sd_rep"]

read_status = snakemake.config["read_status"]
if read_status == 'single' or snakemake.config["seq_tech"] == 'pacbio' or snakemake.config["seq_tech"] == 'ont':
    pair = 1
else:
    pair = 2
if snakemake.config["seq_tech"] == 'illumina':
    rl = snakemake.config["illumina_read_len"]
else:
    rl = snakemake.config["longreads_mean_len"]

"""
Functions
"""


def get_lognormal_dist(tb, mu, sigma):
    d = [random.lognormvariate(mu=mu, sigma=sigma) for i in range(len(tb))]
    tb.insert(loc=tb.shape[1], column='lognormal', value=d)
    tb['RelativeProp'] = tb['lognormal']/tb['lognormal'].sum()
    return tb


def get_even_dist(tb):
    """
    function that calculates even read percentages across superkingdoms
    """
    nb_superkingdoms = len(set(tb['superkingdom']))
    even_prop = 1/nb_superkingdoms
    prop = 0
    props = []
    nb = list(tb['superkingdom']).count('Bacteria')
    na = list(tb['superkingdom']).count('Archaea')
    nv = list(tb['superkingdom']).count('Viruses')
    ne = list(tb['superkingdom']).count('Eukaryota')
    for i in range(len(tb)):
        if tb.iloc[i]['superkingdom'] == 'Bacteria':
            prop = even_prop / nb
        if tb.iloc[i]['superkingdom'] == 'Archaea':
            prop = even_prop / na
        if tb.iloc[i]['superkingdom'] == 'Viruses':
            prop = even_prop / nv
        if tb.iloc[i]['superkingdom'] == 'Eukaryota':
            prop = even_prop / ne
        props.append(prop)
    tb['RelativeProp'] = props
    return tb


def calculate_reads_and_coverage(table, total, sd, read_len, pairing, rep, input_value):
    dic = {}
    genomenames = list(table['AssemblyNames'])
    genomesizes = list(table['Assembly_length'])
    reads_per_genome = 0
    coverage_per_genome = 0
    if input_value == 'RelativeProp' or input_value == 'lognormal' or input_value == 'even':
        props = list(table['RelativeProp'])
        totalgensize = sum(genomesizes)
        fraction_genome = [i/float(totalgensize) for i in genomesizes]
        combined_fractions = [i[0]*i[1] for i in zip(props, fraction_genome)]
        x = total/sum(combined_fractions)
        reads_per_genome = [i*x for i in combined_fractions]
        coverage_per_genome = [(i[0]*read_len*pairing)/i[1] for i in zip(reads_per_genome, genomesizes)]
    if input_value == 'ReadPercent':
        readpercents = list(table['ReadPercent'])
        reads_per_genome = [r*total for r in readpercents]
        coverage_per_genome = [(i[0]*read_len*pairing)/(i[1]) for i in zip(reads_per_genome, genomesizes)]
    if input_value == 'Coverage':
        coverage_per_genome = list(table['Coverage'])
        reads_per_genome = [(i[0]*i[1])/(read_len*pairing)for i in zip(coverage_per_genome, genomesizes)]
    if input_value == 'Reads':
        reads_per_genome = list(table['Reads'])
        coverage_per_genome = [(i[0]*read_len*pairing)/i[1] for i in zip(reads_per_genome, genomesizes)]
    for name, r, c in zip(genomenames, reads_per_genome, coverage_per_genome):
        read = np.random.normal(r, sd)
        cov = np.random.normal(c, sd)
        dic[name] = {}
        dic[name]['AssemblyNames'] = name
        dic[name]['Reads'] = int(read)
        dic[name]['Coverage'] = cov
        logging.info(f"simulating {read} reads, with a coverage value of {cov} reads for accession {name}, "
                     f"replicate {rep} with seed {seed}")
    tb = pd.DataFrame.from_dict(dic, orient='index').reset_index(drop=True)
    return tb


"""
Main
"""
intb = pd.read_csv(snakemake.input.input_tb, sep='\t')
astb = pd.read_csv(snakemake.input.assemblies_tb, sep='\t')
common_col = intb.columns.intersection(astb.columns)[0]
assemblies_with_val = astb.merge(intb, how='left', on=f'{common_col}')
if inputval == 'Coverage' or inputval == 'Reads':
    assemblies_with_val[f'{inputval}'] = np.int64(assemblies_with_val[f'{inputval}']/assemblies_with_val['nb_genomes'])
    assemblies_with_val.drop(columns='nb_genomes', inplace=True)
elif inputval == 'ReadPercent' or inputval == 'RelativeProp':
    assemblies_with_val[f'{inputval}'] = np.float64(assemblies_with_val[f'{inputval}'] /
                                                    assemblies_with_val['nb_genomes'])
    assemblies_with_val.drop(columns='nb_genomes', inplace=True)
elif inputval == 'even':
    assemblies_with_val = get_even_dist(astb)
elif inputval == 'lognormal':
    mu = snakemake.config['mu']
    si = snakemake.config['sigma']
    assemblies_with_val = get_lognormal_dist(astb, mu=mu, sigma=si)
cov_read_tb = calculate_reads_and_coverage(assemblies_with_val, totalreads, sdr, rl, pair, replicate, inputval)
if inputval == 'Coverage' or inputval == 'Reads':
    cov_read_tb = cov_read_tb.drop(inputval, axis=1)
mergedtb = assemblies_with_val.merge(cov_read_tb, how='left', on='AssemblyNames')
rc_tb = mergedtb.drop(['AssemblyInput', 'GbUid', 'FtpPath_GenBank', 'FtpPath_RefSeq', 'AsmReleaseDate_GenBank',
                       'ContigN50', 'ScaffoldN50', 'Assembly_coverage', 'Contig_count'], axis=1)
rc_tb.to_csv(snakemake.output["rc_table"], sep='\t', index=False)
