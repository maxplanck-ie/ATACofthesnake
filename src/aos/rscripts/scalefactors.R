.libPaths(R.home("library"))

suppressMessages(library(edgeR))
suppressMessages(library(tools))

args <- commandArgs(trailingOnly=T)
mat = args[[1]]
outfile = args[[2]]

countmat <- read.csv(mat, sep='\t', header=TRUE)
rowmat <- paste(countmat[,1], countmat[,2], countmat[,3], sep='|')
rownames(countmat) <- rowmat
countmat <- as.matrix(countmat[,-(1:3),drop=FALSE])
# Get norm factors.
NormFactor <- calcNormFactors(object = countmat, method = "TMM")
# Get libSize
LibSize <- colSums(countmat)
# Size factors
SizeFactors <- NormFactor * LibSize / 1000000
# Reciprocal for deepTools
SizeFactors.reci <- 1/SizeFactors
write.table(SizeFactors.reci, file = outfile, col.names = FALSE, quote=FALSE)