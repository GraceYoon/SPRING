# simple example

rm(list = ls())
library(SPRING)

# Load the synthetic count data
data("SynthData2") # n = 1000 and p = 100 synthetic dataset

# SPRING on Synthetic Data, when assuming the data as quantitative counts
\dontrun{
# This takes around 511 seconds
fit.spring <- SPRING(SynthData2, quantitative = TRUE, nlambda = 10, rep.num = 10)
}

# SPRING on Compositional data. Row sums are scaled to 1.
\dontrun{
compoData <- SynthData2/rowSums(SynthData2)
fit.spring <- SPRING(compoData, quantitative = FALSE, nlambda = 10, rep.num = 10)
}
