When using the commands [`mess download`](../commands/download.md) or [`mess run`](../commands/run.md), the pipeline automatically downloads genomes from taxa or accessions in the input table.

## Input 
!!! example
    ```sh
    cat minimal_test.tsv
    taxon	nb	cov_sim	sample
    staphylococcus_aureus	1	0.1	sample1
    1290	1	0.1	sample2
    ```
    See [other examples](input.md#download-table)

## Command

```sh
mess download minimal_test.tsv -o download_test
```

## Output

By default, compressed fasta files are in assembly_finder/download and summary tables are assembly_finder, as highlighted below.

```sh hl_lines="3-7 18 19 21 22"
📂download_test
 ┣ 📂assembly_finder
 ┃ ┣ 📂download
 ┃ ┃ ┣ 📂GCF_000013425.1
 ┃ ┃ ┃ ┗ 📜GCF_000013425.1_ASM1342v1_genomic.fna.gz
 ┃ ┃ ┣ 📂GCF_003812505.1
 ┃ ┃ ┃ ┗ 📜GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┃ ┃ ┗ 📜.snakemake_timestamp
 ┃ ┣ 📂logs
 ┃ ┃ ┣ 📂taxons
 ┃ ┃ ┃ ┣ 📜1290.log
 ┃ ┃ ┃ ┗ 📜staphylococcus_aureus.log
 ┃ ┃ ┣ 📜archive.log
 ┃ ┃ ┣ 📜lineage.log
 ┃ ┃ ┣ 📜rsync.log
 ┃ ┃ ┗ 📜unzip.log
 ┃ ┣ 📜archive.zip
 ┃ ┣ 📜assembly_finder.log
 ┃ ┣ 📜assembly_summary.tsv
 ┃ ┣ 📜config.yaml
 ┃ ┣ 📜sequence_report.tsv
 ┃ ┗ 📜taxonomy.tsv
 ┣ 📂benchmarks
 ┃ ┗ 📂assembly_finder
 ┃ ┃ ┗ 📜download.txt
 ┣ 📂logs
 ┃ ┗ 📜assembly_finder.log
 ┣ 📜config.yaml
 ┣ 📜mess.log
 ┗ 📜samples.tsv
```
