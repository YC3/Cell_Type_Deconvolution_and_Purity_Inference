# Data Purity Estimation (Cell Type Deconvolution)

## Description
In this small project, I present a simple deconvolution method to estimate cell compositions from gene expression profiles of bulk tissues which consists of multiple cell types. 

## Background
Tumor tissues are composed of multiple cell types, including tumor cells and different types of immune cells. Different cell types have different gene expression profiles which makes it difficult to study tumor cell expression patterns.  

Cell type deconvolution methods help us estimate the purity of bulk tumor tissues and also the proportions of cell types in bulk tissues.

Scripts in this small project provide a way to estimate tumor sample heterogeneity by:

1. first modeling the observed mixed gene expression profiles of tumor tissues in linear equations, **_y=Sw_**.
1. and then solving the equations (estimating **_w_**) with quadratic programming (QP).


## A step-by-step guide

### Step 1: Preparing the data

The gene expression (RNA-seq) data of 4 different immune cells were downloaded from Gene Expression Omnibus (GEO) with the accession number 'GSE139242'.

This is a dataset with profiling of human thymic and peripheral blood CD4 + and CD8+ T cells. You can find the description of this dataset [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139242).

Note, this dataset has no tumor cells and is just a mock data used to explain the method. You can directly read it in with: 

```R
data <- read.csv('exp_processed.csv', row.names = 1)
```

### Step 2: Generate a signature matrix for each cell type in the mixture

The cell signature matrix **_S_** is a matrix with rows being cell-type-specific genes, columns being cell types, and entries being averaged gene expression levels in each cell type. **_S_** is usually unknown and need to be estimated from their expression profiles.

Non-cell-type-specific genes show similar expression levels in different types of cell, which contribute little in distinguishing cell types. Thus, this step is to identify cell-type-specific genes that are uniquely expressed in a specific cell type. 

First do pairwise comparisons:

```R
## sample setup (in study GSE139242)
classes <- c(rep('CD4_thymus', 4), rep('CD4_blood', 5), rep('CD8_thymus', 5), rep('CD8_blood', 4))

## pairwise comparisons
comp.o <- pairComp(data, classes, 'fdr')
```

Next, identified top ranked cell-specific genes.

```R
## keep the top 2000 genes in each cell type
n = 2000
CD4_thymus <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD4_thymus', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_blood', n)), selecTop(comp.o$res$'CD4_thymus-CD8_thymus', n))
CD4_blood <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD4_thymus', n), selecTop(comp.o$res$'CD4_blood-CD8_blood', n)), selecTop(comp.o$res$'CD4_blood-CD8_thymus', n))
CD8_thymus <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD8_thymus', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_thymus', n, F)), selecTop(comp.o$res$'CD8_blood-CD8_thymus', n, F))
CD8_blood <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD8_blood', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_blood', n, F)), selecTop(comp.o$res$'CD8_blood-CD8_thymus', n))
```

Then, generate the signature matrix (**_S_**).

```R
feature.v <- unique(c(CD4_thymus, CD4_blood, CD8_thymus, CD8_blood))
source.m <- data[feature.v, ]

S <- t(rmDupRows(t(source.m), classes))
```

### Step 3 Cell type deconvolution 

Create a mock bulk data **_y_** by mixing expression profiles of single cell types in silico.   
You can mix these pure signals in any proportions (sum up to 1) and see how well it estimates your input proportions.

```R
#args <- commandArgs(trailingOnly = TRUE)
args <- c(0, 0.5, 0.4, 0.1)

bulk.v <- S[, "CD4_blood"]*as.numeric(args[1]) + S[, "CD4_thymus"]*as.numeric(args[2]) + S[, "CD8_blood"]*as.numeric(args[3]) + S[, "CD8_thymus"]*as.numeric(args[4]) + rnorm(nrow(S), 0, 0.1)
```

Then, with the mixed bulk data **_y_** and estimated signature matrix **_S_**, we can estimate the proportions of the 4 cell types of cell **_w_** by solving this equation **_y=Sw_** with QP:

```R
## In reality we only has an observed mixed data: bulk.v, and
## an estimated signature data: S.
fraction_infer = purityQP(bulk.v, S)
```

## Want to see results faster?

If you find following the step-by-step guide tedious, trying the command line below would be a good start. 

``` console 
$ Rscript ./PurityInfer.R 0 0.5 0.4 0.1

[1] "estimated cell proportions are:"
[1] 0.0000 0.4999 0.3972 0.1025
```

It does these thing for you:

1. First identifies the significant cell-type-specific genes for each of the 4 cell types;
1. Generates the signature matrix **_S_**;
1. Mix the gene expression profiles of 4 cell types in proportions input by user to mimic the mixed gene expression profile **_y_** in bulk tissues, which is 0%-50%-40%-10% in this example;
1. Estimate the cell proportions of the 'bulk tissue' with QP.

The returned tuple shows the estimated proportions in the 'bulk tissue', which is very close to user input.

