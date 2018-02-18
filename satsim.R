#This function simulates data under a given condition, and runs ML on it, returning the tree and parameter estimates.

satsim <- function(ntaxa = c(16, 32, 64), seqlen = c(500, 1500, 2500), imbal = F, trlen = exp(seq(-3.5, 3.5, 1))[c(1:3, 6:8)], proprandsites = c(0, 0.3, 0.6, 0.9), modelsim = , modelest = ){

}