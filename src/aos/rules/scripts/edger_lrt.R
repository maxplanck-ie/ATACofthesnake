# .libPaths(R.home("library"))
suppressMessages(library("readr"))
suppressMessages(library("edgeR"))

mat = snakemake@input[['mat']]
samplesheet = snakemake@params[['samplesheet']]
of = snakemake@params[['outputfolder']]
edgeR_out = snakemake@output[['table']]
relevantsamples_out = snakemake@output[['samples']]
comparison_entry = snakemake@params[['comparison']]
comparison_name = snakemake@params[['comparison_name']]

get_lrt_coef <- function(full, red, ss) {
    X <- model.matrix(full, data = ss)
    full_terms <- attr(terms(full), "term.labels")
    red_terms  <- attr(terms(red), "term.labels")
    dropped_terms <- setdiff(full_terms, red_terms)
    assign_vec <- attr(X, "assign")
    dropped_idx <- match(dropped_terms, full_terms)
    coef_idx <- which(assign_vec %in% dropped_idx)
    coef_idx
}

# read samplesheet
samplesheet = read.table(
    samplesheet,
    row.names = 1,
    header=TRUE,
    stringsAsFactors=TRUE
)
print("Column classes in samplesheet:")
print(sapply(samplesheet, class))
# read in count matrix.
countmat <- read.csv(mat, sep='\t')
rows <- paste(countmat[,1],countmat[,2], countmat[,3], sep='|')
countmat <- countmat[-c(1:3)]
rownames(countmat) <- rows
# make sure the columns of count matrix are equal to the rows of the samplesheet.
countmat <- countmat[, rownames(samplesheet)] 

# Design
## Either design is present in the comparison_entry, or we default to additive design of all covariates.
if (!is.null(comparison_entry$design)) {
    formula <- as.formula(comparison_entry$design)
    design <- model.matrix(formula, data=samplesheet)
    print(paste("Design for", comparison_name, " given in entry: ", paste(deparse(formula), collapse = "")))
} else {
    formula <- as.formula(
        paste('~',paste(colnames(samplesheet),collapse='+') )
    )
    design <- model.matrix(formula, data=samplesheet)
    print(paste("Design for", comparison_name, " not given in entry. Defaulting to: ", paste(deparse(formula), collapse = "")))
}
## Get coefs for LRT
reduced_formula <- as.formula(comparison_entry$reduced)
lrt_coefs <- get_lrt_coef(full=formula, red=reduced_formula, ss=samplesheet)
lrt_coef_names <- colnames(design)[lrt_coefs]
print(paste("LRT coefs for", comparison_name, ":", paste(lrt_coefs, collapse = ", ")))
print(paste("LRT coef names for", comparison_name, ":", paste(lrt_coef_names, collapse = ", ")))


# Get relevant samples.
mf_full <- model.frame(formula, data = samplesheet, na.action = na.pass)
mf_red  <- model.frame(reduced_formula, data = samplesheet, na.action = na.pass)
ok <- complete.cases(mf_full) & complete.cases(mf_red)
relevantsamples <- rownames(samplesheet)[ok]
write_lines(relevantsamples,relevantsamples_out, sep='\n')

# Run edgeR.
keep <- filterByExpr(
    countmat, design=design,min.count = 5, min.prop = 0.49
)
countmat <- countmat[keep,]
countmat_disp <- estimateGLMCommonDisp(countmat, design, verbose=TRUE)
fit <- glmQLFit(countmat, design=design, dispersion = countmat_disp)
lrt <- glmLRT(fit, coef=lrt_coefs)
res <- data.frame(topTags(lrt, n = Inf))
res$peak_id <- rownames(res)

write.table(
    res,
    file=edgeR_out,
    quote=FALSE,
    sep='\t',
    row.names=FALSE
)