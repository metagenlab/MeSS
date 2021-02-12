import pandas as pd
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
sdfl=snakemake.config["sd_fragment_length"]
reps=snakemake.wildcards.rep

def calculate_reads_from_percent(table,total,sd,rep):
    reads_to_sim={}
    for acc in table.index:
        prop=table.loc[acc]['PercentReads']
        read_nb=int(prop*total)
        read_nb_sample = round(np.random.normal(read_nb, sd))
        # This adds variability for the number of simulated reads between replicates
        reads_to_sim[acc]=read_nb_sample
        logging.info(f"{read_nb_sample} reads for {acc} , replicate {rep}")
    return reads_to_sim

def calculate_reads_from_coverage(table,sd,read_len,pairing,rep):
    reads_to_sim={}
    for acc in table.index:
        coverage=int(table.loc[acc]['coverage'])
        genome_length=int(table.loc[acc]['Assembly_length'])
        read_nb=(coverage*genome_length)/(read_len*pairing)
        read_nb_sample = round(np.random.normal(read_nb, sd))
        # This adds variability for the number of simulated reads between replicates
        reads_to_sim[acc]=read_nb_sample
        logging.info(f"{coverage}X coverage represent {read_nb_sample} reads for {acc} , replicate {rep}")
    return reads_to_sim

"""
Main
"""
table=pd.read_csv(snakemake.input[0],delimiter='\t',index_col='AssemblyNames')
if snakemake.config["use_coverage"]:
    rl=snakemake.config["read_length"]
    read_status = snakemake.config["read_status"]
    if read_status == 'single':
        pair = 1
    else:
        pair = 2
    nb_reads_dic=calculate_reads_from_coverage(table,sdfl,rl,pair,reps)
else:
    total=snakemake.config['total_amount_reads']
    nb_reads_dic=calculate_reads_from_percent(table,total)
table['SimReads']=pd.Series(nb_reads_dic)
table=table.reset_index()
tax_df=table[['AssemblyNames','AssemblyID','AssemblyStatus','Assembly_length','Taxid','superkingdom','phylum','order','family','genus','species','SimReads']]
tax_df_krona=table[['SimReads','superkingdom','phylum','order','family','genus','species']]
tax_df_krona=tax_df_krona.set_index(['SimReads'])#Set all tax ranks as index
tax_df_krona.to_csv(snakemake.output["krona_output"],sep='\t',index=True)
tax_df.to_csv(snakemake.output["rc_table"],sep='\t',index=False)

