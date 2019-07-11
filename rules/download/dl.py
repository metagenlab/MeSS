import os
import ftplib
import pandas as pd



def dl_input(rule_name):
    checkpoint_output=checkpoint.rule_name.get(**wildcards).output[0]
    table=pd.read_csv(checkpoint_output,delimiter='\t')
    acc=list(table['AssemblyAccession'])
    links=list(table['RefseqPath'])


def dl_fna_ftp(ftp_login,full_path,destination):
        print('Downloading to',destination)
        path=full_path.split('ftp://ftp.ncbi.nlm.nih.gov')[1]
        ftp_login.cwd(path)
        filename=snakemake.input[0]
        print('Downloading %s'%filename)
        ftp_login.retrbinary("RETR "+filename, open(os.path.join(destination,filename),"wb").write)
        print('Done')

'''
Main
'''
ftp=ftplib.FTP(host='ftp.ncbi.nih.gov',user='anonymous',passwd=snakemake.params['NCBI_email'])
links=list(pd.read_csv(snakemake.input[0],delimiter='\t')['RefseqPath'])
dl_fna_ftp(ftp,links,snakemake.output[0])
