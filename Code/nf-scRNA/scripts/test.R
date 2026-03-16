args <- commandArgs(trailingOnly=TRUE)
input <- args[1]



writeLines(input, con="results.txt")