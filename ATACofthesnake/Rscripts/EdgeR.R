.libPaths(R.home("library"))

suppressMessages(library("edgeR"))

args <- commandArgs(trailingOnly=T)
mat = args[[1]]
conds = args[[2]]
outfilesign = args[[3]]
outfileall = args[[4]]
batches = args[[5]]

conds = strsplit(conds, ",")[[1]]
countmat <- read.csv(mat, sep='\t')
rows <- paste(countmat[,1],countmat[,2], countmat[,3], sep='_')
countmat <- countmat[-c(1:3)]
rownames(countmat) <- rows
if (length(batches) > 1) {
    design <- model.matrix(~factor(batches, levels(unique(batches))) + factor(conds, levels=(unique(conds))))
} else {
    design <- model.matrix(~factor(conds, levels=(unique(conds))))
}

keep <- filterByExpr(countmat, design=design)
countmat <- countmat[keep,]
countmat_disp <- estimateGLMCommonDisp(countmat, design, verbose=TRUE)
fit <- glmQLFit(countmat, design=design,dispersion = countmat_disp)
res <- glmQLFTest(fit, coef=ncol(fit$design))
ressig <- topTags(res, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
ressig_LFCcut <- ressig$table[abs(ressig$table$logFC) > 1,]
ressig_LFCcut <- ressig_LFCcut[ressig_LFCcut$FDR < 0.05,]
ressig$table$peak_id <- rownames(ressig$table)
ressig_LFCcut$peak_id <- rownames(ressig_LFCcut)
write.table(ressig$table, file=outfileall, quote=FALSE, sep='\t',row.names=FALSE)
write.table(ressig_LFCcut, file=outfilesign, quote=FALSE, sep='\t',row.names=FALSE)