import pandas as pd
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
inputval = snakemake.params.value
if inputval == 'RelativeProp' or inputval == 'ReadPercent':
    totalreads = snakemake.config['total_reads']
else:
    totalreads = 0
sdr = snakemake.config["sd_read_num"]
reps = snakemake.wildcards.rep
read_status = snakemake.config["read_status"]
if read_status == 'single' or snakemake.config["seq_tech"] == 'longreads':
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


def get_even_reads(tb, pv, pb, ph, pe):
    """
    function that calculates the number of viral human and bacterial genomes from the output of assembly_finder
    and assigns even proportions for each genome

    ph=proportion of human reads
    pb=proportion of bacterial reads
    pv=proportion of viral reads
    pe=proportion of eukaryotic reads
    table=output of assembly_finder
    """
    even_prop = []
    nb = list(tb['superkingdom']).count('Bacteria')
    nv = list(tb['superkingdom']).count('Viruses')
    nh = list(tb['species']).count('Homo sapiens')
    ne = list(tb['superkingdom']).count('Eukaryota') - nh  # Eukaryota number, except humans

    for i in range(len(tb)):
        if tb.iloc[i]['superkingdom'] == 'Bacteria':
            even_prop.append(round(pb / nb, 4))#Round to 4 digits after the comma

        if tb.iloc[i]['superkingdom'] == 'Viruses':
            even_prop.append(round(pv / nv, 4))

        if tb.iloc[i]['superkingdom'] == 'Eukaryota' and tb.iloc[i]['species'] == 'Homo sapiens':
            even_prop.append(round(ph / nh, 4))

        if tb.iloc[i]['superkingdom'] == 'Eukaryota' and not tb.iloc[i]['species'] == 'Homo sapiens' and ne > 0:
            even_prop.append(round(pe / ne, 4))

    tb['PercentReads'] = even_prop
    return tb


def calculate_reads_and_coverage(table, total, sd, read_len, pairing, rep, input_value):
    dic = {}
    genomenames = list(table['AssemblyNames'])
    genomesizes = list(table['Assembly_length'])
    reads_per_genome = 0
    coverage_per_genome = 0
    if input_value == 'RelativeProp':
        props = list(table['RelativeProp'])
        totalgensize = sum(genomesizes)
        fraction_genome = [i/float(totalgensize) for i in genomesizes]
        combined_fractions = [i[0]*i[1] for i in zip(props, fraction_genome)]
        x = total/sum(combined_fractions)
        reads_per_genome = [i*x for i in combined_fractions]
        coverage_per_genome = [(i[0]*read_len)/i[1] for i in zip(reads_per_genome, genomesizes)]
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
        dic[name]['Reads'] = read
        dic[name]['Coverage'] = cov
        logging.info(f"simulating {read} reads, with a coverage value of {cov} reads for accession {name}, "
                     f"replicate {rep}")
    tb = pd.DataFrame.from_dict(dic, orient='index').reset_index(drop=True)
    return tb


"""
Main
"""
proportion_reads = {'virus': 0.01, 'human': 0.9, 'bacteria': 0.08, 'non_human_eukaryotes': 0.01}
intb = pd.read_csv(snakemake.input.input_tb, sep='\t')
astb = pd.read_csv(snakemake.input.assemblies_tb, sep='\t')
common_col = intb.columns.intersection(astb.columns)[0]
assemblies_with_val = astb.merge(intb, how='left', on=f'{common_col}')
if inputval == 'Coverage' or inputval == 'Reads':
    assemblies_with_val[f'{inputval}'] = np.int64(assemblies_with_val[f'{inputval}']/assemblies_with_val['nb_genomes'])
    assemblies_with_val.drop(columns='nb_genomes', inplace=True)
elif inputval == 'PercentReads' or inputval == 'RelativeProp':
    assemblies_with_val[f'{inputval}'] = np.float64(assemblies_with_val[f'{inputval}'] /
                                                    assemblies_with_val['nb_genomes'])
    assemblies_with_val.drop(columns='nb_genomes', inplace=True)
else:
    vrp = proportion_reads['virus']
    hrp = proportion_reads['human']
    brp = proportion_reads['bacteria']
    erp = proportion_reads['non_human_eukaryotes']
    assemblies_with_val = get_even_reads(assembly_tb, pv=vrp, pb=brp, ph=hrp, pe=erp)

cov_read_tb = calculate_reads_and_coverage(assemblies_with_val, totalreads, sdr, rl, pair, reps, inputval)
mergedtb = assemblies_with_val.merge(cov_read_tb, on=[f'{inputval}', 'AssemblyNames'])

tax_df = mergedtb[['AssemblyNames', 'RefSeq_category', 'AssemblyStatus', 'Assembly_length', 'Taxid', 'superkingdom',
                   'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'Reads', 'Coverage']]
tax_df.to_csv(snakemake.output["rc_table"], sep='\t', index=False)

