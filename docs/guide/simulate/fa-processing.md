Before generating reads, `mess` pre-processes fastas for read simulators as shown below:

``` mermaid
graph LR
  A[genome_asmV3.fna.gz] --> B{seqkit 
  seq};
  B--> C[genome.fna];
  C --> D{split_contigs.py};
  D --> E[genome_contig1.fna];
  D --> F[genome_contig2.fna];
  D --> G[genome_contig3.fna];
```

## seqkit seq

### File renaming

Files are renamed according to the `fasta` or `accession` column in the input file or assembly_summary.tsv file.

:octicons-arrow-right-24: Done to avoid Snakemake filename matching erros.

### Decompressing

:octicons-arrow-right-24: Done because read simulators take decompressed fasta as input

### Removing long fasta headers

:octicons-arrow-right-24: Done to ensure correct BAM formatting for samtools



## split_contigs.py

### Correct amount of sequences
Art_illumina generates the same amount of reads per contig ([Milhaven and Pfeifer, 2023](https://doi.org/10.1038/s41437-022-00577-3)). This means that, for fragmented assemblies, read sampling will not be uniform across contigs.

In fact, when we use `art_illumina --fcov 10` on this [slightly fragmented E.coli genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000157115.2/), we expect 51534530 simulated bases (10x genome size). However we obtain, 51251100, which less than requested. Thus, to uniformly sample reads acros contigs, we split the genome and apply the same art_illumina command for each contig. 

### Compatibilty with pbsim3
Since pbsim3 already splits fasta, inputting split fasta is faster.

