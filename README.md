
<!-- README.md is generated from README.Rmd. Please edit that file -->
SPRING
======

<!-- badges: start -->
<!-- badges: end -->
The R package `SPRING` (Semi-Parametric Rank-based approach for INference in Graphical model) estimates sparse microbial association networks using rank-based correlation with sparse graphical modeling techniques. The corresponding reference is

Yoon, Gaynanova and Müller (2019) [Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data](https://www.frontiersin.org/articles/10.3389/fgene.2019.00516/full). *Frontiers in Genetics*, 10:516.

Installation
------------

``` r
# install.packages("devtools")
devtools::install_github("GraceYoon/SPRING")
```

Example
-------

``` r
library(SPRING)
data("SynthData2") # load saved synthetic count data in this package, of dimension n = 1000 and p = 100

# SPRING on Synthetic Data, when assuming the data as quantitative counts
# This takes around 511 seconds
fit.spring <- SPRING(SynthData2, quantitative = TRUE, nlambda = 10, rep.num = 10)

# SPRING on Compositional data. Row sums are scaled to 1.
compoData <- SynthData2/rowSums(SynthData2)
fit.spring <- SPRING(compoData, quantitative = FALSE, nlambda = 10, rep.num = 10)
```
