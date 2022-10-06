.libPaths(R.home("library"))

suppressMessages(library("edgeR"))

args <- commandArgs(trailingOnly=T)
mat = args[[1]]
conds = args[[2]]
outfilesign = args[[3]]
outfileall = args[[4]]
batches = args[[5]]

conds = unlist(strsplit(conds, ","))
countmat <- read.csv(mat, sep='\t')
rows <- paste(countmat[,1],countmat[,2], countmat[,3], sep='_')
countmat <- countmat[-c(1:3)]
rownames(countmat) <- rows

condF =  factor(conds, levels=(unique(conds)))

batches = unlist(strsplit(batches, ","))
print(conds)
print(batches)
if (length(unique(batches)) > 1) {
    batchF = factor(batches, levels=(unique(batches)))
    design <- model.matrix(~batchF + condF)
} else {
    design <- model.matrix(~condF)
}

keep <- filterByExpr(countmat, design=design,min.count = 5, min.prop = 0.49)
countmat <- countmat[keep,]
countmat_disp <- estimateGLMCommonDisp(countmat, design, verbose=TRUE)
fit <- glmQLFit(countmat, design=design,dispersion = countmat_disp)
res <- glmQLFTest(fit, coef=ncol(fit$design))
ressig <- topTags(res, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
#ressig_LFCcut <- ressig$table[abs(ressig$table$logFC) > 1,]
ressig_LFCcut <- ressig$table[ressig$table$FDR < 0.05,]
ressig$table$peak_id <- rownames(ressig$table)
ressig_LFCcut$peak_id <- rownames(ressig_LFCcut)
write.table(ressig$table, file=outfileall, quote=FALSE, sep='\t',row.names=FALSE)
write.table(ressig_LFCcut, file=outfilesign, quote=FALSE, sep='\t',row.names=FALSE)