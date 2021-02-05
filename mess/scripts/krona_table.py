import pandas as pd
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)

def calculate_reads(table,total):
    reads_to_sim=[]
    for i in range(len(table)):
        acc=table['AssemblyNames'][i]
        prop=table.loc[i]['PercentReads']
        read_nb=int(prop*total)
        mean = read_nb
        sd = snakemake.config["sd_read_num"]
        read_nb_sample = round(np.random.normal(mean, sd))
        # This adds variability for the number of simulated reads between replicates
        reads_to_sim.append(read_nb_sample)
        logging.info(f"{read_nb_sample} reads for {acc} , replicate {snakemake.wildcards.rep}")
    return reads_to_sim
"""
Main
"""
checkpoint_output = snakemake.input[0]
table=pd.read_csv(checkpoint_output,delimiter='\t')
total=snakemake.config['total_amount_reads']
nb_reads_list=calculate_reads(table,total)
tax_df_index=table.loc[:,['Taxid','superkingdom','phylum','order','family','genus','species']]
tax_df_no_index=table.loc[:,['superkingdom','phylum','order','family','genus','species']]
tax_df_index.index=table['AssemblyNames']
tax_df_index.insert(loc=0,column='SimReads',value=nb_reads_list)
tax_df_no_index.insert(loc=0,column='SimReads',value=nb_reads_list)
tax_df_no_index.to_csv(snakemake.output["no_index"],sep='\t',index=False)
tax_df_index.to_csv(snakemake.output["index"],sep='\t',index=True)

