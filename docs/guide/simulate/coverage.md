Before generating reads, [`mess simulate`](../commands/simulate.md) calculates coverage depths for each genome according to values in the input table or a distribution set by the user.

If replicates are set, and by default, the same coverage values will be applied for each replicate. This can be modified by increasing the standard deviation between each replicates with `--rep-sd`.

## From Input table

``` mermaid
graph LR
  A(total bases);
  B[bases] -->| genome size| D[coverage depth];
  C[reads] --> |read length x pairing| B[bases];
  D[coverage depth] --> F[simulator];
```


## From distributions



