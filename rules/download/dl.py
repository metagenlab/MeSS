import os
import ftplib
import pandas as pd

def dl_fna_ftp(ftp_login,wildcards):
        checkpoint_output = checkpoints.generate_table_summary.get(**wildcards).output[0]
        table=pd.read_csv(checkpoint_output,delimiter='\t',index_col=0)
        full_path=table[wildcards.sample]['RefseqPath']
        path=full_path.split('ftp://ftp.ncbi.nlm.nih.gov')[1]
        ftp_login.cwd(path)
        filename=wildcards.sample
        print('Downloading %s'%filename)
        ftp_login.retrbinary("RETR "+filename, open(os.path(snakemake.output[0]),"wb").write)
        print('Done')

'''
Main
'''
ftp=ftplib.FTP(host='ftp.ncbi.nih.gov',user='anonymous',passwd=snakemake.params['NCBI_email'])

dl_fna_ftp(ftp,wildcards)
