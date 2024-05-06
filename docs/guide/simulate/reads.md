`mess` feeds each [processed genome](fa-processing.md) with its corresponding [coverage depth](coverage.md) to [art_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) or [pbsim3](https://github.com/yukiteruono/pbsim3) to generate fastq files and alignment files (optionally).


``` mermaid
flowchart TB
  fa[fasta]
  fq[fastq]
  sam[sam]
  bam[bam]
  bai[bam.bai]
  maf[maf]
  cov[coverage.txt]
  tax[taxonomic 
  profile]
  fqgz[fastq.gz]
  
  art(art_illumina)
  cat(cat)
  pigz(pigz)
  pbsim(pbsim3)
  bioconvert{bioconvert}
  sam2bam(sam2bam)
  maf2sam(maf2sam)
  seqkit{seqkit}
  shufle(shuffle)
  replace(replace)
  samtools{samtools}
  sortidenx(sort & index)
  coverage(coverage)
  taxonkit(taxonkit 
  profile2cami)

  fa --> art
  fa --> pbsim
  art --> fq
  pbsim --> fq
  art --> sam
  pbsim --> maf

  maf --> bioconvert
  sam --> bioconvert
  sam2bam --> bam
  bam --> samtools
  samtools --> sortindex
  sortindex --> bai
  samtools --> coverage
  coverage --> cov
  cov --> taxonkit
  taxonkit --> tax
  bioconvert --> sam2bam
  bioconvert --> maf2sam
  
  fq --> cat
  cat --> pigz
  pigz --> fqgz
```