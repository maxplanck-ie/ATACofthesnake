.libPaths(R.home("library"))

suppressMessages(library(edgeR))
suppressMessages(library(tools))

mat = snakemake@input[[1]]
outfile = snakemake@output[[1]]

countmat <- read.csv(mat, sep='\t', header=TRUE)
rowmat <- paste(countmat[,1], countmat[,2], countmat[,3], sep='|')
rownames(countmat) <- rowmat
countmat <- as.matrix(countmat[,-(1:3),drop=FALSE])
# Get norm factors.
NormFactor <- calcNormFactors(object = countmat, method = "TMM")
# Get libSize
LibSize <- colSums(countmat)
if (any(LibSize == 0)) {
  bad <- colnames(countmat)[LibSize == 0]
  stop(sprintf(
    "Sample(s) %s have zero reads overlapping the peak set; cannot compute a scale factor for them. Check upstream filtering/peak calling for these samples.",
    paste(bad, collapse = ", ")
  ))
}
# Size factors
SizeFactors <- NormFactor * LibSize / 1000000
# Reciprocal for deepTools
SizeFactors.reci <- 1/SizeFactors
write.table(SizeFactors.reci, file = outfile, col.names = FALSE, quote=FALSE)