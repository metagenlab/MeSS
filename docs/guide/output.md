## Directory structure
```sh hl_lines="3-7 15-24 28-30"
📂mess_out
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

## Coverage table

Check genome size (total_sequence_length), coverage depth (cov_sim) and other stats for both samples in cov.tsv

=== "sample1"
    | fasta | total_sequence_length | number_of_contigs | reads | bases | seq_abundance | cov_sim | tax_abundance |
    | --------------- | --------------------- | ----------------- | ------- | -------- | ------------- | ------- | ------------- |
    | GCF_000013425.1 | 2821361 | 1 | 940.454 | 282136.1 | 1.0 | 0.1 | 1.0 |

=== "sample2"
    | fasta           | total_sequence_length | number_of_contigs | reads   | bases    | seq_abundance | cov_sim | tax_abundance |
    | --------------- | --------------------- | ----------------- | ------- | -------- | ------------- | ------- | ------------- |
    | GCF_003812505.1 | 2257431               | 3                 | 752.477 | 225743.1 | 1.0           | 0.1     | 1.0           |

## Reads

Let's check the fastq content using [`seqkit stats`](https://bioinf.shenwei.me/seqkit/usage/#stats).

```sh
$ seqkit stats --all test_bam_out/fastq/*.fq.gz
file                    format  type  num_seqs  sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)
fastq/sample1_R1.fq.gz  FASTQ   DNA        940  141,000      150      150      150  150  150  150        0  150        1   97.99   91.72    27.68  32.89
fastq/sample1_R2.fq.gz  FASTQ   DNA        940  141,000      150      150      150  150  150  150        0  150        1   97.35    89.7    26.52  33.19
fastq/sample2_R1.fq.gz  FASTQ   DNA        752  112,800      150      150      150  150  150  150        0  150        1   98.02   91.75    27.78  31.38
fastq/sample2_R2.fq.gz  FASTQ   DNA        752  112,800      150      150      150  150  150  150        0  150        1   97.27    89.7    26.51  31.52
```

- [x] 90% of reads are Q30
- [x] Requested amount of sequences were simulated (940 and 752 reads for sample1 and sample2 respectively)
- [x] Both genomes are covered at 0.1x (divide sum_len by genome size)


## Alignments

Using the `--bam` flag, sorted and indexed BAM files are generated, which can be quality controlled with [qualimap](http://qualimap.conesalab.org/). Qualimap reports for each sample can be aggregated using [multiQC ](https://multiqc.info/).

### Qualimap
```sh
cd mess_out/bam 
for bam in *.bam; do qualimap bamqc -bam $bam; done
```
### MultiQC
```sh
multiqc .
```

##### General stats


| Sample        | QualiMap_mqc-generalstats-qualimap-avg_gc | QualiMap_mqc-generalstats-qualimap-median_insert_size | QualiMap_mqc-generalstats-qualimap-1_x_pc | QualiMap_mqc-generalstats-qualimap-5_x_pc | QualiMap_mqc-generalstats-qualimap-10_x_pc | QualiMap_mqc-generalstats-qualimap-30_x_pc | QualiMap_mqc-generalstats-qualimap-50_x_pc | QualiMap_mqc-generalstats-qualimap-median_coverage | QualiMap_mqc-generalstats-qualimap-mean_coverage | QualiMap_mqc-generalstats-qualimap-general_error_rate | QualiMap_mqc-generalstats-qualimap-percentage_aligned | QualiMap_mqc-generalstats-qualimap-mapped_reads | QualiMap_mqc-generalstats-qualimap-total_reads |
| ------------- | ----------------------------------------- | ----------------------------------------------------- | ----------------------------------------- | ----------------------------------------- | ------------------------------------------ | ------------------------------------------ | ------------------------------------------ | -------------------------------------------------- | ------------------------------------------------ | ----------------------------------------------------- | ----------------------------------------------------- | ----------------------------------------------- | ---------------------------------------------- |
| sample1       |                                           |                                                       |                                           |                                           |                                            |                                            |                                            |                                                    | 0.1                                              | 0.0                                                   | 100.0                                                 | 1880                                            | 1880                                           |
| sample1_stats | 33.48098434004474                         | 199                                                   | 6.449724087062946                         | 0.0                                       | 0.0                                        | 0.0                                        | 0.0                                        | 0                                                  |                                                  |                                                       |                                                       |                                                 |                                                |
| sample2       |                                           |                                                       |                                           |                                           |                                            |                                            |                                            |                                                    | 0.0999                                           | 0.0                                                   | 100.0                                                 | 1504                                            | 1504                                           |
| sample2_stats | 31.822164948453608                        | 199                                                   | 6.446044198028644                         | 0.0003543851395679425                     | 0.0                                        | 0.0                                        | 0.0                                        | 0                                                  |                                                  |                                                       |                                                       |                                                 |                                                |

##### qualimap bamqc genome stats

| Sample  | total_reads | mapped_reads | mapped_bases | sequenced_bases | duplication_rate | mean_insert_size | median_insert_size | mean_mapping_quality | gc_percentage | general_error_rate | homopolymer_indels | mean_coverage | percentage_aligned |
| ------- | ----------- | ------------ | ------------ | --------------- | ---------------- | ---------------- | ------------------ | -------------------- | ------------- | ------------------ | ------------------ | ------------- | ------------------ |
| sample1 | 1880        | 1880         | 281999       | 281998          | 0.05             | 199.8862         | 199.0              | 91.0798              | 32.96         | 0.0                | 0.0                | 0.1           | 100.0              |
| sample2 | 1504        | 1504         | 225601       | 225600          | 0.0              | 199.3763         | 199.0              | 83.9774              | 31.39         | 0.0                | 0.0                | 0.0999        | 100.0              |

- [x] 0% general error rate and homopolymer indels
- [x] At least 83 mapping quality with 100% read alignment
- [x] Average coverage depth values are of 0.1x

!!! warning
    Only 6% of both genome locations are covered at 1x. This is to be expected as our requested coverage depth is at 0.1x. If we want to cover a larger fraction of the genomes, we need to simulate more reads.

## Taxonomic profiles

=== "sample1"
    ```sh
    @SampleID:sample1
    @Version:0.10.0
    @Ranks:superkingdom|phylum|class|order|family|genus|species|strain
    @TaxonomyID:
    @@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
    2	superkingdom	2	Bacteria	100.000000000000000
    1239	phylum	2|1239	Bacteria|Bacillota	100.000000000000000
    91061	class	2|1239|91061	Bacteria|Bacillota|Bacilli	100.000000000000000
    1385	order	2|1239|91061|1385	Bacteria|Bacillota|Bacilli|Bacillales	100.000000000000000
    90964	family	2|1239|91061|1385|90964	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae	100.000000000000000
    1279	genus	2|1239|91061|1385|90964|1279	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus	100.000000000000000
    1280	species	2|1239|91061|1385|90964|1279|1280	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus|Staphylococcus aureus	100.000000000000000
    93061	strain	2|1239|91061|1385|90964|1279|1280|93061	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus|Staphylococcus aureus|Staphylococcus aureus subsp. aureus NCTC 8325	100.000000000000000
    ```

=== "sample2"
    ```sh
    @SampleID:sample2
    @Version:0.10.0
    @Ranks:superkingdom|phylum|class|order|family|genus|species|strain
    @TaxonomyID:
    @@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
    2	superkingdom	2	Bacteria	100.000000000000000
    1239	phylum	2|1239	Bacteria|Bacillota	100.000000000000000
    91061	class	2|1239|91061	Bacteria|Bacillota|Bacilli	100.000000000000000
    1385	order	2|1239|91061|1385	Bacteria|Bacillota|Bacilli|Bacillales	100.000000000000000
    90964	family	2|1239|91061|1385|90964	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae	100.000000000000000
    1279	genus	2|1239|91061|1385|90964|1279	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus	100.000000000000000
    1290	species	2|1239|91061|1385|90964|1279|1290	Bacteria|Bacillota|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus|Staphylococcus hominis	100.000000000000000
    ```
- [x] Samples have the correct taxonomic profile in biobox format