# Input types

All mess commands require as input a tsv file or a directory containing tsv files, specified with `--input`.

## Directory

The sequencing_run directory contains 3 tsv files, one for each metagenomic community.

```sh
📂sequencing_run
 ┣ 📜sample1.tsv
 ┣ 📜sample2.tsv
 ┗ 📜sample3.tsv
```

## TSV
There are three types of columns in the input table:

* taxon/accession
* nb
* abundance (reads, bases, tax_abundance, seq_abundancea and cov_sim)


The combination of taxon/accession and nb columns, is used by [`mess download`](../commands/download.md) to download a number of taxa or accessions.

The abundance column is used by [`mess simulate`](../commands/simulate.md) to calculate the coverage depth of each genome.

### Download table 

!!! example 

    === "taxon"

        | taxon                  | nb  |
        | :--------------------- | :-- |
        | 1280  | 1   |
        | pseudomonas_aeruginosa | 1   |

    === "accession"

        | accession              |
        | :--------------------- |
        | GCF_000418345.1 |
        | GCF_000233495.1 |

### Simulation table

!!! example 
    [`mess simulate`](../commands/simulate.md) uses coverage depth as input for read simulators.

    Reads, bases, taxonomic or sequence abundance will be [converted to coverages](community.md).

    === "coverage"

        | taxon                  | nb  | cov_sim |
        | :--------------------- | :-- | :------ |
        | 1280  | 1   | 10      |
        | pseudomonas_aeruginosa | 1   | 10      |

    === "bases"

        | taxon                  | nb  | bases    |
        | :--------------------- | :-- | :------- |
        | 1280  | 1   | 28213610 |
        | pseudomonas_aeruginosa | 1   | 62644040 |

    === "reads"

        | taxon                  | nb  | reads  |
        | :--------------------- | :-- | :----- |
        | 1280  | 1   | 94045  |
        | pseudomonas_aeruginosa | 1   | 208813 |

    === "tax_abundance"

        | taxon                  | nb  | tax_abundance  |
        | :--------------------- | :-- | :----- |
        | 1280  | 1   | 0.5  |
        | pseudomonas_aeruginosa | 1   | 0.5 |

    === "seq_abundance"

        | taxon                  | nb  | seq_abundance  |
        | :--------------------- | :-- | :----- |
        | 1280  | 1   | 0.32  |
        | pseudomonas_aeruginosa | 1   | 0.68 |

### Local genomes table

!!! example 
    [`mess simulate`](../commands/simulate.md) needs genome path, size and tax_ids in the input table.

    | fasta  | path                   | cov_sim | tax_id | total_sequence_length |
    | :----- | :--------------------- | :------ | :----- | :-------------------- |
    | fasta1 | /path/to/fasta1.fna.gz | 10      | taxid1      | 6666666          |
    | fasta2 | /path/to/fasta2.fna.gz | 10      | taxid2      | 5555555          |
