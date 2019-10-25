
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
data("QMP") # load the data available from this package.

# Apply SPRING on QMP data.
fit.spring <- SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific", nlambda = 50, rep.num = 50)
# This takes around 23 minutes. We are working on reducing the computation time (10/25/2019).

# StARS-selected lambda index based on the threshold (default = 0.01)
opt.K <- fit.spring$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K <- fit.spring$fit$est$path[[opt.K]]
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K <- SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode = 'maxabs')
```
