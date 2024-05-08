`mess` feeds [split fastas](fa-processing.md) and their [coverage depths](coverage.md) to [art_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) or [pbsim3](https://github.com/yukiteruono/pbsim3) for read simulation as shown below.


``` mermaid
flowchart TB
  %%file types
  fa[fasta]
  fq[fastq]
  fqgz[fastq.gz]
  sam[sam]
  maf[maf]
  bam[bam]
  bai[bam <br/>bam.bai]
  cov[coverage.txt]
  tax[biobox<br/>taxonomic<br/>profile]

  %%tools
  art(art_illumina)
  pbsim(pbsim3)
  cat(cat)
  pigz(pigz)
  sam2bam(bioconvert <br/>sam2bam)
  maf2sam(bioconvert <br/>maf2sam)
  shuffle(seqkit <br/>shuffle)
  replace(seqkit <br/>replace)
  merge(samtools <br/>merge)
  sortindex(samtools <br/>sort & index)
  coverage(samtools <br/>coverage)
  taxonkit(taxonkit <br/>profile2cami)

  %%workflow 
  fa --> art
  fa --> pbsim
  fa2tax[genome taxids]
  subgraph simulators
    art
    pbsim
  end 
  art --> fq
  art ==>|contigs sam| sam
  pbsim --> fq
  pbsim ==>|no header contigs <br/>maf| sed(sed)
  sed ==>|add contig header| maf
  pbsim -.-> ccs

  fq -->|contigs<br/>fastq| pigz
  pigz -->|contigs<br/>fastq.gz| cat
  cat -->|samples fastq| shuffle
  shuffle -->|shuffled fastq| replace
  replace -->|anonymized fastq| fqgz

  maf ==> maf2sam
  maf2sam ==> sam
  sam ==> sam2bam
  sam2bam ==> bam
  bam ==>|contigs bam| merge
  merge ==>|samples bam| sortindex
  sortindex ==>|sorted bam| coverage
  fa2tax ==> cov
  coverage ==>|contigs<br/>coverages| cov
  cov ==>|taxids<br/>coverages| taxonkit

  ccs[ccs.sam] -.-> sam2bam
  bam -.-> pbccs(pbccs)
  pbccs -.->|hifi <br/>fastq.gz<br/>contigs| cat

  taxonkit ==> tax
  sortindex ==>|sorted<br/>indexed<br/>bam| bai  
  
  subgraph output
    fqgz
    bai
    tax 
  end
```

>:material-arrow-right: Path for fastq files

>:material-arrow-right-bold: Path for alignments (sam, bam, maf...)

>**---**> Path for PacBio hifi fastq

## Reads processing
Simulators output uncompressed reads for every contig, which are then compressed with pigz and concatenated with cat.
Concatenated fastq are then shuffled and anonymized with seqkit shuffle and seqkit replace respectively. 

## Alignments processing
Simulators can optionally output alignment files, usually in SAM format, which are then converted to sorted and indexed BAM.
In some case, alignment files are in other formats and need some processing to be correctly conterted to BAM. For example, pbsim3 output MAF alignment files with no sequence ID. `mess` which will add sequence IDs, and convert MAF to SAM and finally BAM using bioconvert. 