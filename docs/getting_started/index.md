# Getting started
## Quick start 
### [Install with mamba](../install.md#mamba-recommended)

```sh
mamba create -n mess mess
```
### Run 
```sh
mess test --bam
```

### Output
```sh hl_lines="3-7 15-24 28-30"
📂test_out
 ┣ 📂assembly_finder
 ┃ ┣ 📂download
 ┃ ┃ ┣ 📂GCF_000418345.1
 ┃ ┃ ┃ ┗ 📜GCF_000418345.1_ASM41834v1_genomic.fna.gz
 ┃ ┃ ┣ 📂GCF_003812505.1
 ┃ ┃ ┃ ┗ 📜GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┃ ┃ ┗ 📜.snakemake_timestamp
 ┃ ┣ 📜archive.zip
 ┃ ┣ 📜assembly_finder.log
 ┃ ┣ 📜assembly_summary.tsv
 ┃ ┣ 📜config.yaml
 ┃ ┣ 📜sequence_report.tsv
 ┃ ┗ 📜taxonomy.tsv
 ┣ 📂bam
 ┃ ┣ 📜sample1.bam
 ┃ ┣ 📜sample1.bam.bai
 ┃ ┣ 📜sample2.bam
 ┃ ┗ 📜sample2.bam.bai
 ┣ 📂fastq
 ┃ ┣ 📜sample1_R1.fq.gz
 ┃ ┣ 📜sample1_R2.fq.gz
 ┃ ┣ 📜sample2_R1.fq.gz
 ┃ ┗ 📜sample2_R2.fq.gz
 ┣ 📂processing
 ┃ ┣ 📜cov.tsv
 ┃ ┗ 📜split.tsv
 ┣ 📂tax
 ┃ ┣ 📜sample1_profile.txt
 ┃ ┗ 📜sample2_profile.txt
 ┣ 📜config.yaml
 ┗ 📜mess.log
```

## [Step-by-step guide](../guide/index.md)