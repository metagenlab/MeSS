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


tb_w_reads=get_even_reads(snakemake.input[0],pv=vrp,pb=brp,ph=hrp,pe=erp)

tb_w_reads.to_csv(snakemake.output[0], sep='\t', header=True)