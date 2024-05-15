# PBSIM3
You can customize `pbsim3` parameters to generate a wide range of Nanopore error profiles

### R10.4.1

By default, `mess` uses the high quality error profile (QSHMM-ONT-HQ) based on the R10.4 flowcell, which simulates reads [with similar accuracy compared to the R10.4.1 flowcell](https://github.com/yukiteruono/pbsim3/issues/12).

In addition, `mess` uses the pbsim3 recommended Nanopore substitution:insertion:deletion ratio (39:24:36) coupled with 0.99 accuracy as show below:


```yaml
ratio: '39:24:36'
max_len: 1000000
mean_len: 9000
min_len: 100
sd_len: 7000
accuracy: 0.99
model: QSHMM-ONT-HQ
```

Simulate R10.4.1 reads with the test command

```sh
mess test \
    --tech nanopore \
    --output r10.4.1
```

### R10.3

You can tweak `mess` parameters to mimic older Nanopore flowcells

```sh
mess test \
    --tech nanopore \
    --model QSHMM-ONT \
    --accuracy 0.95 \
    --output r10.3
```