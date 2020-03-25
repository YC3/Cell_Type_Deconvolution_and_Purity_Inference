# Data Purity Estimation

Functions in this repository can be used to estimate data purity and proportions of subtypes.

The gene expression data of blood cells were downloaded from Gene Expression Omnibus which was used as pure signals. You could mix these pure signals in any proportions (sum up to 1) and see how well it estimates your input.

## Example
```
$ Rscript ./PurityInfer.R 0 0.5 0.4 0.1

[1] "estimated cell proportions are:"
[1] 0.0000 0.4999 0.3972 0.1025
```
