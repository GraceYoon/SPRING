## code to prepare `DATASET` dataset goes here

# To generate synthetic data using extended american gut data.

rm(list=ls())

library(SpiecEasi) # make_graph and graph2prec function.
library(mixedCCA) # calculate Kendall Correlation. function "estimateR"
source("../copulaMicrobiome/Rfunctions/synthData.R") # this has the function "synthData_from_ecdf"

seed = 10010
# load the amgut real data (subset of real data: size is n=2000 and p=p1 with minimum depth 1e4)
load(paste0("../copulaMicrobiome/Data/amgutsim_p", p1, ".rdata"))


n <- 500; p1 <- 200
    set.seed(seed)
    gtype = "scale_free"
    e1 = 2*p1 # number of edges

    set.seed(seed) # set the seed number for make_graph part.
    graph_p1 <- SpiecEasi::make_graph(gtype, p1, e1)
    Prec1  <- SpiecEasi::graph2prec(graph_p1)
    Cor1   <- cov2cor(prec2cov(Prec1))

    # True counts
    SynthData <- synthData_from_ecdf(get(paste0("amgutsim_p", p1)), Sigma = Cor1, n = n, seed = seed)

    usethis::use_data(SynthData)


n <- 1000; p1 <- 100
    set.seed(seed)
    gtype = "cluster"
    e1 = 2*p1 # number of edges

    set.seed(seed) # set the seed number for make_graph part.
    graph_p1 <- SpiecEasi::make_graph(gtype, p1, e1)
    Prec1  <- SpiecEasi::graph2prec(graph_p1)
    Cor1   <- cov2cor(prec2cov(Prec1))

    # True counts
    SynthData2 <- synthData_from_ecdf(get(paste0("amgutsim_p", p1)), Sigma = Cor1, n = n, seed = seed)

    usethis::use_data(SynthData2)




load("../copulaMicrobiome/Analysis/SimpleEx/qmphealthyrank6pruned.rdata")
# this containing three data variable:
# X: copyadjusted count data
# QMP: quantitative count data (X -> RMP by dividing by total abundance -> QMP by multiplying by cell counts)
# qmphealthy6_only1filt: phyloseq class data.

### But only save QMP for this package.
usethis::use_data(QMP)


