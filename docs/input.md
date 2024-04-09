MeSS requires an input tsv file or a directory containing tsv files, specified with `--input`.

## Directory

```
📦gut_subsets
 ┣ 📜sample1.tsv
 ┣ 📜sample2.tsv
 ┗ 📜sample3.tsv
```

## TSV

The tsv file requires an accession/taxon column to download genomes and a value to calculate the coverage for metagenome simulation.

### Download table

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

You can provide reads, bases, taxonomic or sequence abundance which will be converted to coverages.

Alternatively, you can provide the coverage value directly.

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

To simulate from local fasta, MeSS needs the genomes tax_ids and sequence lengths in the input table.

| fasta  | path                   | cov_sim | tax_id | total_sequence_length |
| :----- | :--------------------- | :------ | :----- | :-------------------- |
| fasta1 | /path/to/fasta1.fna.gz | 10      | taxid1      | 6666666          |
| fasta2 | /path/to/fasta2.fna.gz | 10      | taxid2      | 5555555          |
