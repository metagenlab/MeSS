# PBSIM3
You can customize `pbsim3` parameters to generate PacBio CLR or HIFI reads

### CLR
Defaults:

```yaml
ratio: '22:45:33'
max_len: 1000000
mean_len: 9000
min_len: 100
sd_len: 7000
accuracy: 0.85
model: QSHMM-RSII
```

Simulate CLR reads with the test command

```sh
mess test \
    --tech pacbio \
    --output pacbio_clr
```

### HIFI

You can add `--error hifi` to simulate highly accurate PacBio reads with 10 passes

Defaults:
```yaml
ratio: '22:45:33'
max_len: 1000000
mean_len: 9000
min_len: 100
sd_len: 7000
accuracy: 0.999
passes: 10
model: QSHMM-RSII
```

Run:
```sh
mess test \
    --tech pacbio \
    --error hifi \
    --output pacbio_hifi
```