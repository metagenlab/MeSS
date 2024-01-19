# `MeSS` Input

## Input

MeSS requires an input tsv file or directory, specified with `--input`.

### Directory

```
📦gut_subsets
 ┣ 📜sample1.tsv
 ┣ 📜sample2.tsv
 ┗ 📜sample3.tsv
```

### Input table columns

the tsv file has to have an entry and nb headers.

To set each entry's quantity, you can provide a reads, bases, coverage or abundance column. Optionally, the table can have a sample column, if not, the sample name will be the file's basename.

Examples:

- Coverage

| entry                  | nb  | cov_sim |
| :--------------------- | :-- | :------ |
| staphylococcus_aureus  | 1   | 10      |
| pseudomonas_aeruginosa | 1   | 10      |

- Bases

| entry                  | nb  | bases    |
| :--------------------- | :-- | :------- |
| staphylococcus_aureus  | 1   | 28213610 |
| pseudomonas_aeruginosa | 1   | 62644040 |

- Reads

| entry                  | nb  | reads  |
| :--------------------- | :-- | :----- |
| staphylococcus_aureus  | 1   | 94045  |
| pseudomonas_aeruginosa | 1   | 208813 |
