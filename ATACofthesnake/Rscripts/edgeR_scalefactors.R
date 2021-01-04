.libPaths(R.home("library"))

library(edgeR)
library(tools)
args <- commandArgs(trailingOnly=T)
mat = args[[1]]

countmat <- read.csv(mat, sep='\t', header=TRUE)
head(countmat)
touch('diffAcc/scalefactors.txt')