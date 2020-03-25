#!/usr/local/bin Rscript
args = commandArgs(trailingOnly = TRUE)

library(limma)
library(quadprog)

# Author: YC3
#
# You could change input proportions to mix signal from 4 sources, and compare with the output result.
#
# Example: Rscript ./PurityInfer.R 0.25 0.25 0.3 0.2



pairComp <- function(data, classes, padj){
    
    # design matrix
    classFac.v = as.factor(classes)
    classLev.v <- levels(classFac.v)
    design.m <- model.matrix(~ 0 + classFac.v)
    colnames(design.m) <- classLev.v

    # pair-wise contrast matrix
    N <- length(classLev.v)
    compsNum <- 0.5*N*(N - 1)

    colN.v <- vector()
    for(i in 1:(N - 1)) colN.v <- c(colN.v, paste0(classLev.v[i], '-', classLev.v[(i + 1):N]))
    contr.m <- matrix(0, nrow = N, ncol = compsNum, dimnames = list(classLev.v, colN.v))

    for (i in 1:(N - 1)){

        clab <- classLev.v[(i + 1):N] 
        colN.v0 <- paste0(classLev.v[i], '-', clab)

        contr.m[classLev.v[i], colN.v0] <- 1
        for (j in clab) contr.m[j, colN.v0[clab == j]] <- -1
    }
    
    # linear model fit and empirical Bayesian estimation 
    # of adjusted variance and differential expression
    fit.o <- contrasts.fit(lmFit(data, design.m), contr.m)
    ebayes.o <- eBayes(fit.o)

    # select top genes for every conparison
    res.l <- vector(length = compsNum, mode = 'list')
    names(res.l) <- colnames(contr.m)
    for(i in 1:compsNum) res.l[[i]] <- topTable(ebayes.o, coef = i, adjust.method = padj, number = nrow(data))

    list(res = res.l, contrast = contr.m)
}



selecTop <- function(resTable, topN, selPos = TRUE){
    
    genes.v = rownames(resTable)
    if(selPos) sig.v <- genes.v[intersect(which(resTable[, 't'] > 0), which(resTable[, 'adj.P.Val'] < 0.05))]
    else sig.v <- genes.v[intersect(which(resTable[, 't'] < 0), which(resTable[, 'adj.P.Val'] < 0.05))]

    if(length(sig.v) < topN) sig.v
    else sig.v[1:topN]
}



purityQP <- function(x, S){

    D.m <- 2 * t(S) %*% S
    d.v <- 2 * x %*% S # or ignore 2 in both D and d.
    A.m <- diag(ncol(S))
    b.v <- rep(0, ncol(S)) # cell proportions for each cell type >= 0.

    qp.o <- solve.QP(D.m, d.v, A.m, b.v)
    qp.o$sol

}



rmDupRows <- function(data, classFac){
    sum.o <- summary(as.factor(classFac), maxsum = length(unique(classFac)))
    sum.m <- rowsum(data, group = classFac)
    data.new <- sweep(sum.m, 1, sum.o, "/")
  
    data.new
}



# read in data
data = read.csv('exp_processed.csv', row.names = 1)

# cell-type specific gene
classes <- c(rep('CD4_thymus', 4), rep('CD4_blood', 5), rep('CD8_thymus', 5), rep('CD8_blood', 4))
comp.o <- pairComp(data, classes, 'fdr')

# merge samples according to PCA results
n = 2000
CD4_thymus <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD4_thymus', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_blood', n)), selecTop(comp.o$res$'CD4_thymus-CD8_thymus', n))
CD4_blood <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD4_thymus', n), selecTop(comp.o$res$'CD4_blood-CD8_blood', n)), selecTop(comp.o$res$'CD4_blood-CD8_thymus', n))
CD8_thymus <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD8_thymus', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_thymus', n, F)), selecTop(comp.o$res$'CD8_blood-CD8_thymus', n, F))
CD8_blood <- intersect(intersect(selecTop(comp.o$res$'CD4_blood-CD8_blood', n, F), selecTop(comp.o$res$'CD4_thymus-CD8_blood', n, F)), selecTop(comp.o$res$'CD8_blood-CD8_thymus', n))

feature.v <- unique(c(CD4_thymus, CD4_blood, CD8_thymus, CD8_blood))
source.m <- data[feature.v, ]
S <- t(rmDupRows(t(source.m), classes))

# Generating mixed data

bulk.v <- S[, "CD4_blood"]*as.numeric(args[1]) + S[, "CD4_thymus"]*as.numeric(args[2]) + S[, "CD8_blood"]*as.numeric(args[3]) + S[, "CD8_thymus"]*as.numeric(args[4]) + rnorm(nrow(S), 0, 0.1)

# Estimating cell proportions in bulk tissue
print('estimated cell proportions are:')
print(round(purityQP(bulk.v, S), 4))
