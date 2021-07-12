import pandas as pd
"""
Main
"""
ftplinks = pd.read_csv(snakemake.input[0], sep='\t')['FtpPath_GenBank']
links = []
for link in ftplinks:
    link = link.replace('ftp://ftp.ncbi.nlm.nih.gov', '')
    fna = '/' + link.split('/')[-1]+'_genomic.fna.gz\n'
    link += fna
    links.append(link)

f = open(snakemake.output[0], "w")
f.writelines(links)
f.close()
