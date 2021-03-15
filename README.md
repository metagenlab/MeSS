#MeSS [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/version.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/latest_release_date.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/downloads.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/platforms.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/license.svg)](https://anaconda.org/metagenlab/mess)
The Metagenomic Sequence Simulator (MeSS) is a snakemake workflow used for simulating metagenomic mock communities.
## Installation
[![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/installer/conda.svg)](https://conda.anaconda.org/metagenlab)    
```shh
conda install -c conda-forge mamba
mamba install -c metagenlab -c conda-forge -c bioconda mess
```
## Required files
### Input table example
MeSS takes the same input as Assembly_finder, with an additional column for either Coverage values, read percentages or
relative abundances.
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
community_name: metagenome-1
proportion_reads: {'virus':0.01,'human':0.9,'bacteria':0.08,'non_human_eukaryotes':0.01}

#Sequencing run params
seq_tech: longreads #[illumina,longreads]
read_status: single 
total_reads: 100000

##Illumina (art params)
illumina_sequencing_system: HSXt #HiSeqX TruSeq (read length:150bp)
illumina_read_len: 150
illumina_mean_frag_len: 200
illumina_sd_frag_len: 10

##Long reads (pbsim2 params)
chemistry: R94 
longreads_min_len: 100
longreads_max_len: 1000000
longreads_sd_len: 7000
longreads_mean_len: 9000
longreads_mean_acc: 85
difference_ratio: "23:31:46"

## Replicate params
replicates: 1 
sd_read_num: 0 #Standard deviation for varying the number of reads per replicate, if 0 keep the same number of reads for each replicate
set_seed: 1 #Set this parameter to a fixed number so data from random number generators is reproducible

#Assembly finder parameters
input_table_path: input_table.tsv
NCBI_key: your_ncbi_key
NCBI_email: your_ncbi_email
complete_assemblies: False
reference_assemblies: False
representative_assemblies: False
exclude_from_metagenomes: True
Genbank_assemblies: True
Refseq_assemblies: True
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

## Running MeSS

Here is an example command to run MeSS on the previously described config and input table.
```bash
snakemake --snakefile path/to/MeSS/Snakefile --configfile config.yml --use-conda --conda-prefix path/to/conda/envs/ --resources ncbi_requests=3 nb_simulation=2 parallel_cat=2 --cores 10 illumina_sim
```

**nb_simulation** controlls the number of parallel art_illumina jobs. This parameter is important for memory usage, as art_illumina loads genomes into memory. Thus, for big genomes it is recommended to lower this parameter.

**parallel_cat** controlls the number of genomes to be concatenated in parallel. For big genomes and computers with low memory, lowering this parameter lowers memory usage.

## MeSS outputs
### Directory structure
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

