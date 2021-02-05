import pandas as pd
def get_even_reads(table, pv, pb, ph, pe):
    '''
    function that calculates the number of viral human and bacterial genomes from the output of assembly_finder
    and assigns even proportions for each genome

    ph=proportion of human reads
    pb=proportion of bacterial reads
    pv=proportion of viral reads
    pe=proportion of eukaryotic reads
    table=output of assembly_finder
    '''
    even_prop = []
    tb = pd.read_csv(table, delimiter='\t', index_col=0)
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

vrp = snakemake.params['proportion_reads']['virus']
hrp = snakemake.params['proportion_reads']['human']
brp = snakemake.params['proportion_reads']['bacteria']
erp = snakemake.params['proportion_reads']['non_human_eukaryotes']

tb = pd.read_csv(snakemake.input.table, sep='\t')
# Check if the user already inputted read percentages for each assembly
if 'PercentReads' in list(tb.columns) and len(list(tb['PercentReads'])) == len(tb):
    # Check if the user filled the list or not
    Percent_reads_list = True
    assembly_table=pd.read_csv(snakemake.params['table'],sep='\t')
    tb_w_reads = pd.merge(assembly_table,tb,how='outer',on='UserInputNames')
    tb_w_reads['PercentReads'] = tb_w_reads['PercentReads'] / tb_w_reads['nb_genomes']#If multiple assemblies per entry, divide read percent by the number of genomes
    tb_w_reads.drop('nb_genomes', axis=1, inplace=True)
    tb_w_reads.set_index('AssemblyNames',inplace=True)

else:#If no read percent is specified
    Percent_reads_list = False
    tb_w_reads=get_even_reads(snakemake.params['table'],pv=vrp,pb=brp,ph=hrp,pe=erp)

tb_w_reads.to_csv(snakemake.output[0], sep='\t', header=True)
