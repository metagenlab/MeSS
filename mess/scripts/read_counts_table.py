import pandas as pd
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
inputval=snakemake.params.value
totalreads = snakemake.config['total_reads']
sdr=snakemake.config["sd_read_num"]
reps=snakemake.wildcards.rep
read_status = snakemake.config["read_status"]
if read_status == 'single' or snakemake.config["seq_tech"]=='longreads':
    pair = 1
else:
    pair = 2
if snakemake.config["seq_tech"]=='illumina':
    rl=snakemake.config["illumina_read_len"]
else:
    rl=snakemake.config["longreads_mean_len"]


def calculate_reads_and_coverage(table,total,sd,read_len,pairing,rep,input_value):
    dic={}
    genomenames=list(table['AssemblyNames'])
    genomesizes=list(table['Assembly_length'])
    if input_value=='RelativeProp':
        props=list(table['RelativeProp'])
        totalgensize=sum(genomesizes)
        fraction_genome = [i/float(totalgensize) for i in genomesizes]
        combined_fractions = [i[0]*i[1] for i in zip(props, fraction_genome)]
        x = total/sum(combined_fractions)
        reads_per_genome = [i*x for i in combined_fractions]
        coverage_per_genome = [(i[0]*read_len)/i[1] for i in zip(reads_per_genome, genomesizes)]
    if input_value=='ReadPercent':
        readpercents=list(table['ReadPercent'])
        reads_per_genome=[r*total for r in readpercents]
        coverage_per_genome=[(i[0]*read_len*pairing)/(i[1]) for i in zip(reads_per_genome,genomesizes)]
    if input_value=='Coverage':
        coverage_per_genome=list(table['Coverage'])
        reads_per_genome=[(i[0]*i[1])/(read_len*pairing)for i in zip(coverage_per_genome,genomesizes)]
    for name,r,c in zip(genomenames,reads_per_genome,coverage_per_genome):
        read=round(np.random.normal(r, sd))
        cov=round(np.random.normal(c, sd))
        dic[name]={}
        dic[name]['SimReads']=read
        dic[name]['Coverage']=cov
        logging.info(f"simulating {read} reads, with a coverage value of {cov} reads for accession {name}, replicate {rep}")
    tb=pd.DataFrame.from_dict(dic,orient='index')
    return tb

"""
Main
"""
table=pd.read_csv(snakemake.input[0],sep='\t')
cov_read_tb=calculate_reads_and_coverage(table,totalreads,sdr,rl,pair,reps,inputval)
table.set_index('AssemblyNames',inplace=True)
mergedtb=table.join(cov_read_tb).reset_index()
tax_df=mergedtb[['AssemblyNames','AssemblyID','AssemblyStatus','Assembly_length','Taxid','superkingdom','phylum','order','family','genus','species','SimReads','Coverage',f'{inputval}']]
krona=mergedtb[['SimReads','superkingdom','phylum','order','family','genus','species']].set_index('SimReads')
krona.to_csv(snakemake.output["krona_output"],sep='\t',index=True)
tax_df.to_csv(snakemake.output["rc_table"],sep='\t',index=False)

