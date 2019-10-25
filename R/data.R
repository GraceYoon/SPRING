#' Synthetic count data
#' @name SynthData
#' @description  SynthData and SynthData2 were generated using empirical cdf of American Gut Project Data. SynthData has scale_free-type-graph structure of size 500 rows and 200 columns, and SynthData2 has cluster-type-graph structure of size 1000 rows and 100 columns.
#'
#' @source
#'
#' Yoon, Gaynanova and Müller (2019) Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data. \emph{Frontiers in Genetics.} 10:516. \url{doi:10.3389/fgene.2019.00516}
#' @format \code{SynthData} is an object of class \code{matrix} with 500 rows and 200 columns. \code{SynthData2} is an object of class \code{matrix} with 1000 rows and 100 columns.
"SynthData"



#' @name SynthData
#' @aliases SynthData2
"SynthData2"

#' Quantitative Microbiome Project data
#'
#' @description  The data containing quantitative microbiome count data of dimenstion 106 samples/subjects (in rows) and 91 OTUs (in columns). The raw dataset is pruned the taxa present less than 30\% of samples and final dataset contains only healthy subjects from two cohorts: Study cohort and Disease cohort.
#'
#' @source
#'
#' Yoon, Gaynanova and Müller (2019) Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data. \emph{Frontiers in Genetics}. 10:516. \url{doi:10.3389/fgene.2019.00516}
#'
#' Vanderputte et al. (2017) Quantitative microbiome profiling links gut community variation to microbial load. \emph{Nature}. 551: 507-511. \url{doi:10.1038/nature24460}
#'
"QMP"
