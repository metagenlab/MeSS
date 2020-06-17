
# Metagenomic Sequence Simulator
MeSS is a snakemake workflow used for simulating metagenomic mock communities.
## Dependencies
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.8.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
## Installation
#### Conda
Download and install miniconda 3 (Linux 64-bit)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Snakemake
Install mamba and use it to create the snakemake environment as the conda solver is slow and [struggles to select the latest version](https://github.com/conda/conda/issues/9905)
```bash
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
``` 


#### Clone repository
MeSS contains assembly finder as a submodule, paste the command below to clone both repos:
```bash
git clone --recursive https://github.com/metagenlab/MeSS.git
```

## Required files
### Input table example
MeSS takes the same input as Assembly_finder, with an additional column for the percentage of total simulated read per query.
```
UserInputNames           nb_genomes     PercentReads     
1813735                  1              0.3
114185                   1              0.4
ATCC_13985               3              0.3
```
If the percent read column is not present, MeSS will generate an even distribution within superkingdoms.
In the input table shown above, all queries are bacteria and each query will have a read percentage of 20% (1/5).
### Config file example
```
#MeSS parameters
replicates: 1
community_name: 'sample-metagenome'
proportion_reads: {'virus':0,'human':0,'bacteria':1,'non_human_eukaryotes':0}
##art_illumina parameters
illumina_sequencing_system: 'HSXt'#HiSeqX TruSeq 
total_amount_reads: 10000 #reads to simulate in total
read_length: 150
mean_fragment_length: 200
sd_fragment_length: 10
read_status: 'single'
sd_read_num: 0 #Standard deviation for varying the number of reads per replicate

#Assembly_finder parameters
##E-utilities parameters
input_table_path: path/to/input_table
NCBI_key: your_ncbi_api_key
NCBI_email: your_ncbi_email
##Search filter parameters
complete_assemblies: False
reference_assemblies: False
representative_assemblies: False
exclude_from_metagenomes: True
Genbank_assemblies: True
Refseq_assemblies: True
##Filtering function parameter
Rank_to_filter_by: 'None'
```
#### MeSS parameters
MeSS offers the possibility to generate multiple mock communities using the same set of assembly files in the same directory.
For this, the user has to set up one configuration file per mock community and change the community_name accordingly.
Also, a set of replicates can be generated per one community by changing the replicates parameter.
The user can also control the read percentage attributed to each superkingdom by changing percentage values in the proportion_reads dictionary.
#### art_illumina parameters

MeSS uses [art_illumina](https://academic.oup.com/bioinformatics/article/28/4/593/213322) to generate reads, and the user can change parameters like read length and fragment length under the art_illumina params section as shown above. 
#### Assembly_finder parameters
MeSS uses [Assembly_finder](https://github.com/metagenlab/assembly_finder) to download genomes, and requires the user to have an NCBI account. For more details on Assembly_finder parameters check its documentation.

##Running MeSS
Here is an example command to run MeSS on the previously described config and input table.
```bash
snakemake --snakefile path/to/MeSS/Snakefile --configfile config.yml --use-conda --conda-prefix path/to/conda/envs/ --resources ncbi_requests=3 nb_simulation=2 parallel_cat=2 --cores 10 all_sim 
```

**nb_simulation** controlls the number of parallel art_illumina jobs. This parameter is important for memory usage, as art_illumina loads genomes into memory. Thus, for big genomes it is recommended to lower this parameter.

**parallel_cat** controlls the number of genomes to be concatenated in parallel. For big genomes and computers with low memory, lowering this parameter lowers memory usage.

##MeSS outputs
###Directory structure
After running MeSS, your working directory should look like this:
```bash
├── assembly_gz
│   ├── assembly-accession-1.fna.gz
│   └── assembly-accession-2.fna.gz
├── krona
│   ├── sample-1-metagenome
│   │   └── replicate-1-krona-plot.html 
│   └── sample-2-metagenome
│       └── replicate-1-krona-plot.html
├── logs
│   ├── downloads
│   ├── filtered
│   ├── not-filtered  
│   ├── read_counts_table
│   ├── read_generation
│   └── shuffling
├── simreads
│   ├── sample-1-metagenome
│   │    └── replicate-1
│   │        └── simulated-reads.fastq
│   ├── sample-2-metagenome
│   │    └── replicate-1
│   │        └── simulated-reads.fastq
├── tables
│   ├── filtered
│   ├── not-filtered  
│   └── table-for-krona
├── config.yaml
├── input_table.tsv
└── metagenome-summary.tsv
```
The simulated reads fastqs are compressed and located in the 'simreads' directory, and their composition is summarized in 
'metagenome-summary.tsv'. 

