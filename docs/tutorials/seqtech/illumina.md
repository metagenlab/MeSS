# ART
`art_illumina` is the default illumina read simulator in `mess` and can generate reads with built-in or custom error profiles.

## Built-in error profile

from `art_illumina` help:
```sh
GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)
```
By default, `mess` uses the Hiseq 2500 error model to generate paired 150 bp read lengths with a mean and standard deviation fragment length of 200 and 10 respectively.

Built-in error profile, mean read length, mean and standard deviation fragment length can be changed with `--error`, `--mean-length`, `--frag-len` and `--frag-sd` respectively.


### HiSeq 2000

```sh
mess test \
    --tech illumina \
    --error HS20 \
    --mean-len 100 \
    --frag-len 200 \
    --frag-sd 10 \
    --output hs20
```

### Miseq v3

```sh
mess test \
    --tech illumina \
    --error MSv3 \
    --mean-len 250 \
    --frag-len 500 \
    --frag-sd 10 \
    --output msv3
```

## Custom error profiles
If you want to use updated error profiles, your best bet is to use custom profiles, which can be generated using [art_illumina_profiles](https://github.com/hzi-bifo/art_illumina_profiles). 

!!! tip
    `art_illumina_profiles` needs alignments from reference data like [Genome in a bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle). You can find bam and bam.bai in [NCBI's reference samples ftp ](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/)



### GIAB hiseq 300x

Download forward and reverse error profiles
```sh
curl -O https://raw.githubusercontent.com/hzi-bifo/art_illumina_profiles/master/example_profile/RMNISTHS_30xdownsample.1.txt \
    -O https://raw.githubusercontent.com/hzi-bifo/art_illumina_profiles/master/example_profile/RMNISTHS_30xdownsample.2.txt
```
Run test command with custom profile (maximum read length is 148bp)

```sh
mess test \
    --tech illumina \
    --custom-err RMNISTHS_30xdownsample \
    --mean-len 148 \
    --frag-len 300 \
    --frag-sd 10 \
    --output giab_hs300x
```

### MBARC-26
`CAMISIM` also provides custom pre-built error profiles like the [MBARC-26 error profile](https://www.nature.com/articles/sdata201681), which can be used by `mess`.

Download forward and reverse error profiles
```sh
curl -O https://raw.githubusercontent.com/CAMI-challenge/CAMISIM/master/tools/art_illumina-2.3.6/profiles/ART_MBARC-26_HiSeq_R1.txt \
    -O https://raw.githubusercontent.com/CAMI-challenge/CAMISIM/master/tools/art_illumina-2.3.6/profiles/ART_MBARC-26_HiSeq_R2.txt
```

Run test command with custom profile

```sh
mess test \
    --tech illumina \
    --custom-err ART_MBARC-26_HiSeq_R \
    --mean-len 150 \
    --frag-len 270 \
    --frag-sd 27 \
    --output mbarc26
```