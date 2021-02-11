.libPaths(R.home("library"))

suppressMessages(library("DESeq2"))
suppressMessages(library(tools))

args <- commandArgs(trailingOnly=T)
mat = args[[1]]
conds = args[[2]]
outfilesign = args[[3]]
outMAplot= args[[4]]
outfileall = args[[5]]



conds = strsplit(conds, ",")[[1]]

countmat <- read.csv(mat, sep='\t')
rows <- paste(countmat[,1],countmat[,2], countmat[,3], sep='_')
countmat <- countmat[-c(1:3)]
rownames(countmat) <- rows
cols <- cbind(colnames(countmat), conds)
colnames(cols) <- c("Sample", "Cond")
cols <- as.data.frame(cols)
matkeep = apply(countmat, 1, function(x) all(x > 1))
countmat = countmat[matkeep,]
rows <- rownames(countmat)
dds <- DESeqDataSetFromMatrix(countData = countmat,
                                   colData = cols,
                                   design = ~ Cond)
dds = DESeq(dds)

res <- as.data.frame(results(dds, pAdjustMethod = "fdr" ))
res <- res[complete.cases(res), ]
ressig <- res[res$padj < 0.01,]
write.csv(ressig, file=outfilesign, quote=FALSE, sep='\t')
write.csv(res, file=outfileall, quote=FALSE, sep='\t')