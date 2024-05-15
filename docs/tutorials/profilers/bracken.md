## Input

### Download an example

We can use a bracken report [from bracken's sample data](https://github.com/jenniferlu717/Bracken/blob/master/sample_data/sample_output_species_abundance.txt).

```sh
curl https://raw.githubusercontent.com/jenniferlu717/Bracken/master/sample_data/sample_output_species_abundance.txt -o bracken.txt
```

### Format table

We can use [csvtk](https://github.com/shenwei356/csvtk) to extract the required columns (taxonomy_id and new_est_reads) for the `mess` input table.

```sh
csvtk cut -t -f taxonomy_id,new_est_reads bracken.txt \
    | csvtk rename -t -f 1,2 -n taxon,reads \
    | csvtk -t sort -k reads:nr > bracken.tsv
```

### Update taxids

Use taxonkit to update removed/merged taxids

```sh
csvtk cut -t -U -f taxon bracken.tsv \
    | taxonkit lineage --show-status-code 2> tax.log \
    | csvtk -t add-header -n old_taxid,updated_taxid,lineage 1> lineage.tsv
```

Checking the tax.log (as of may 2024)

```sh
$ cat tax.log
11:03:42.688 [WARN] taxid 1172 was merged into 264691
11:03:42.688 [WARN] taxid 86662 was merged into 1405
```

Get the updated taxids, with taxid count

```sh
csvtk join -t bracken.tsv lineage.tsv -f 1 \
    | csvtk cut -t -f updated_taxid,reads \
    | csvtk rename -t -f 1,2 -n taxon,reads > updated.tsv
```

### Replace unclassified taxids

!!! failure
    Some taxids have no genome data available !

    For example, 60481, 164757 and 189918 correspond to unclassified Shewanella sp. MR-7, Mycobacterium sp. JLS and Mycobacterium sp. KMS, repsectively. 
    
    We can replace unclassified species with other species from the same genus. 

    We also add the nb column to tell mess to download one genome per taxon.

```sh
cat updated.tsv \
    | csvtk replace -t -f taxon -p "60481" -r "211586" \
    | csvtk replace -t -f taxon -p "164757" -r "83332" \
    | csvtk replace -t -f taxon -p "189918" -r "557599" \
    > clean.tsv
```

### Subsample

For a quick run (~ 5 min with 10 threads), we can subset the top 10 most abundant taxa for our simulation.

```sh
csvtk head -t -n 10 clean.tsv > subsample.tsv
```

??? example "subsample.tsv"
    | taxon  | reads   |
    | ------ | ------- |
    | 562    | 1449813 |
    | 1396   | 1288824 |
    | 1280   | 1212241 |
    | 1076   | 824296  |
    | 83332  | 699935  |
    | 148447 | 656706  |
    | 132919 | 622314  |
    | 272131 | 590648  |
    | 303    | 515978  |
    | 32046  | 466493  |

## Simulate

Now that we cleaned up and subsampled bracken's taxonomic profile we can use it for `mess run`.

By default, mess downloads 1 reference genome per taxon and excludes atypical genomes (see [download guide](../guide/download.md).
)
For read simulation, the default parameters generate paired 150 bp illumina reads using art_illumina's HiSeq 2500 error profile.

Using the `--bam` flag yields an alignment file and a ground truth taxonomic profile. 
```sh
mess run -i subsample.tsv --bam --threads 10
```

## Output
### fastq
```sh
$ seqkit stats mess_out/fastq/*.fq.gz --all
processed files:  2 / 2 [======================================] ETA: 0s. done
file                               format  type   num_seqs        sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  Q20(%)  Q30(%)  AvgQual  GC(%)
mess_out/fastq/subsample_R1.fq.gz  FASTQ   DNA   8,326,992  1,249,048,800      150      150      150  150  150  150        0  150   98.01   91.67    27.79  50.76
mess_out/fastq/subsample_R2.fq.gz  FASTQ   DNA   8,326,992  1,249,048,800      150      150      150  150  150  150        0  150   97.31   89.64    26.51  50.76
```

- [x] Total read count correspond to the one requested in the input (256 reads difference) 

### Taxonomic profile
Bracken composition

| taxon  | reads   | seq_abundance        |
| ------ | ------- | -------------------- |
| 562    | 1449813 | 0.17410469821482438  |
| 1396   | 1288824 | 0.15477190063271803  |
| 1280   | 1212241 | 0.14557522485219607  |
| 1076   | 824296  | 0.09898780485461704  |
| 83332  | 699935  | 0.08405357928573762  |
| 148447 | 656706  | 0.07886230841209485  |
| 132919 | 622314  | 0.07473225248005103  |
| 272131 | 590648  | 0.07092955559868037  |
| 303    | 515978  | 0.06196260757455525  |
| 32046  | 466493  | 0.056020068094525345 |


Simulated composition 

I calculated sequence abundances from true read counts in `logs/coverage/subsample.tsv`, and added each genomes taxid at species level (with `taxonkit reformat --show-lineage-taxids`).

| taxon  | reads   | seq_abundance      |
| ------ | ------- | ------------------ |
| 562    | 1449799 | 0.1741083695048584 |
| 1396   | 1288786 | 0.1547720953736955 |
| 1280   | 1212160 | 0.1455699729265982 |
| 1076   | 824275  | 0.0989883261566721 |
| 1773   | 699930  | 0.0840555629211604 |
| 148447 | 656686  | 0.0788623310794582 |
| 132919 | 622291  | 0.0747317879013213 |
| 272131 | 590622  | 0.0709286138379861 |
| 303    | 515970  | 0.0619635517843658 |
| 32046  | 466473  | 0.0560193885138835 |

I used a wilcoxon test to compare abundances
```python
>>> from scipy.stats import wilcoxon
>>> delta = abs(real_df["seq_abundance"] - sim_df["seq_abundance"])
>>> wilcoxon(delta)
WilcoxonResult(statistic=24.0, pvalue=0.76953125)
```

- [x] Non significant difference in sequence abundance between the simulated and real sample (mean delta = ${1.47}^{-06}$, $p=0.77$) 