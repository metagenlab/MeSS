# MeSS [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/version.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/latest_release_date.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/downloads.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/platforms.svg)](https://anaconda.org/metagenlab/mess) [![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/license.svg)](https://anaconda.org/metagenlab/mess)
The Metagenomic Sequence Simulator (MeSS) is a snakemake workflow used for simulating metagenomic mock communities.
## Installation
```bash
git clone --recursive https://github.com/metagenlab/MeSS.git
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
### Or

[![Anaconda-Server Badge](https://anaconda.org/metagenlab/mess/badges/installer/conda.svg)](https://conda.anaconda.org/metagenlab)    
```bash
conda install -c conda-forge mamba
mamba install -c metagenlab -c conda-forge -c bioconda mess
```
## Required files
### Input table examples
MeSS takes the same input as Assembly_finder, with an additional column for either coverage values, read percentages or
relative abundances.
#### Read percentage
Below is an example of input table where the user can set, for each entry, read percentages of the total metagenomic reads

UserInputNames | nb_genomes | PercentReads
--- | --- | ---
1813735 | 1 | 0.3
114185 | 1 |  0.4
ATCC_13985 | 3 |  0.3

If the percent read column is not present, MeSS will generate an even distribution within superkingdoms.
In the input table shown above, if no PercentReads is present, each entry will have a read percentage of 20%
(as all entries belong to the same superkingdom: bacteria)
#### Coverage values
The user has also the option to set coverage values instead of %reads of the total metagenomic reads for each entry.

UserInputNames | nb_genomes | Coverage
--- | --- | ---
1813735 | 1 | 20
114185 | 1 |  30
ATCC_13985 | 3 |  20

In this case, all 3 assemblies found for ATCC_13985 will have the same coverage value of 20
#### Relative proportions 
Alternatively, the user can specify relative proportions between assemblies. Given the total number of reads to
be present in the metagenome, scripts will calculate coverage and read numbers respecting the relative proportions.

UserInputNames | nb_genomes | RelativeProp
--- | --- | ---
1813735 | 1 | 0.3
114185 | 1 |  0.4
ATCC_13985 | 3 |  0.3

For ATCC_13985, the 3 genomes will have a RelativeProp value of 0.1.
### Config file example
```yaml
#MeSS parameters
input_table_path: input_table.tsv
community_name: metagenome-1
#Replicates parameters
replicates: 1 
sd_read_num: 0 

#Random seeds
seed: 1 

#Sequencing run params
seq_tech: ont #[illumina, ont, illumina]
read_status: single 
total_reads: 100000

#Illumina (art params)
illumina_sequencing_system: HSXt #HiSeqX TruSeq (read length:150bp)
illumina_read_len: 150
illumina_mean_frag_len: 200
illumina_sd_frag_len: 10

#Long reads (pbsim2 params)
chemistry: R94 
longreads_min_len: 100
longreads_max_len: 1000000
longreads_sd_len: 7000
longreads_mean_len: 9000
longreads_mean_acc: 85
difference_ratio: "23:31:46"

#Assembly download
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
#### Mess parameters
The path to the input table can be set by the input_table_path parameter in the config file as shown above.

MeSS offers the possibility to generate multiple mock communities using the same set of assembly files in 
the same directory.
For this, the user has to set up one configuration file per mock community and change the community_name accordingly. 
#### Replicates parameters
The user has the option to te create a set of replicates for one community. Each replicate read number can be drawn 
from a normal distribution with a standard deviation set in the sd_read_num parameter. 
#### Random seeds
The MeSS workflow uses random seeds for read generation and read shuffling. To ensure reproducible results, one can
give the seed parameter a fixed number.

#### Sequencing run params
MeSS offers the possibility to select art_illumina or pbsim2 to simulate illumina and long reads respecitvely.
In addition, read pairing and the total amount of reads can be set using the read_status and total_reads parameters.

#### Illumina (art params)
MeSS uses [art_illumina](https://academic.oup.com/bioinformatics/article/28/4/593/213322) to generate illumina reads, 
and the user can change parameters like read and fragment length under the art_illumina params section as shown the 
yaml file above. 

#### Long reads (pbsim2 params)
For long read simulation, [pbsim2](https://github.com/yukiteruono/pbsim2) was integrated in the pipeline.
pbsim2 randomly samples reads from a reference sequence following a gamma distribution, and errors are introduced following 
FIC-HMM models for different chemistries.

For a PC64 PacBio simulation, the user can choose the PC64.model and a difference_ratio of 6:50:54
(substitution:insertion:deletion) as recommended by pbsim2's help message.
As for a Nanopore sequencing run using a R9.4 flowcell, the user can set the values shown in the yaml file above.

For more details check [pbsim2's documentation](https://github.com/yukiteruono/pbsim2/blob/master/README.md)

#### Assembly download
MeSS uses [Assembly_finder](https://github.com/metagenlab/assembly_finder) to download genomes, and requires the user 
to have an NCBI account. For more details on Assembly_finder parameters check its documentation.

## Running MeSS
### Snakemake command
Here is an example command to run MeSS on the previously described config and input table.
```bash
snakemake --snakefile path/to/MeSS/Snakefile --configfile config.yml \
--use-conda --conda-prefix path/to/conda/envs/ \
--resources ncbi_requests=3 nb_simulation=2 parallel_cat=2 --cores 10 all_sim
```

**nb_simulation** controlls the number of parallel art_illuminaor pbsim2 jobs. 
This parameter is important for memory usage, as art_illumina loads genomes into memory.
Thus, for big genomes it is recommended to lower this parameter.

**parallel_cat** controlls the number of genomes to be concatenated in parallel. 
For big genomes and computers with low memory, lowering this parameter lowers memory usage.
### Using the wrapper
```bash
mess run -c 10 
```
Runs the Mess workflow using 10 cores.
## MeSS outputs
### Directory structure
After running MeSS, your working directory should look like this:
```
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
│   └── not-filtered  
├── config.yaml
├── input_table.tsv
├── summary-metagenome-1.tsv
├── summary-metagenome-2.tsv
└── assemblies-summary.tsv
```
The simulated reads fastqs are compressed and located in the simreads/ directory, and their composition is summarized in 
summary-<community_name>.tsv. 

