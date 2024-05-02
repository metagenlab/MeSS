
All [mess commands](../commands/index.md) take as input one or multiple tables.

In this guide we are using [`mess test`](index.md#command) to execute all the workflow steps using a minimal example (shown below)
``` title="minimal_test.tsv"
--8<-- "mess/test_data/minimal_test.tsv"
```


The input file contains:

* Taxon/accession and nb columns for genome download 
* Third column to calculate genome coverage for read simulation.
* sample column to map genomes to their respective sample



## Other examples
### Table
Instead of coverage depths, you can set bases, reads, taxonomic or sequence abundance in the input table. 

See [coverage calculation](simulate/coverage.md) for more details
!!! example
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
    
    === "coverage"

        | taxon                  | nb  | cov_sim |
        | :--------------------- | :-- | :------ |
        | 1280  | 1   | 10      |
        | pseudomonas_aeruginosa | 1   | 10      |


### Directory

If you want to run the pipeline on multiple samples you can point to a directory with multiple tables (one for each sample, with the sample name in the file name).
!!! example 
    ```sh
    📂sequencing_run
    ┣ 📜sample1.tsv
    ┣ 📜sample2.tsv
    ┗ 📜sample3.tsv
    ```
Or you can aggregate all sample info in one table
 
!!! example
    | taxon                 | nb  | cov_sim | sample  |
    | :-------------------- | :-- | :------ | :------ |
    | staphylococcus_aureus | 1   | 10     | sample1 |
    | 1290                  | 1   | 10     | sample2 |
    | 562                  | 1   | 10     | sample3 |





