from Bio import SeqIO

try :
    record=SeqIO.read(snakemake.input[0],'fasta')#Check if the fasta file has one record only
    SeqIO.write(record,snakemake.output[0],'fasta')#If only one record (i.e one contig) copy the file
except ValueError:#More than one record, we merge contigs
    fasta=open(snakemake.input[0],'r')
    handle=SeqIO.parse(fasta,'fasta')
    headers=[]
    seq=[]
    for record in handle:
        headers.append(record.description)
        seq.append(record.seq)
    nlist='N'*1000#Seperate contigs with a string of 1000 N
    merged_contigs=[seq[i]+nlist for i in range(len(seq))]
    merged_contigs[0]='>'+headers[0]+'\n'+merged_contigs[0]#add the first header before the first contig
    with open(snakemake.output[0],'w') as f:
        for contig in merged_contigs:
            f.write(f"{contig}\n")