# .libPaths(R.home("library"))
suppressMessages(library("readr"))
suppressMessages(library("edgeR"))

mat = snakemake@input[['mat']]
samplesheet = snakemake@params[['samplesheet']]
pseudocount = as.integer(snakemake@params[['pseudocount']])
of = snakemake@params[['outputfolder']]
edgeR_out = snakemake@output[['table']]
relevantsamples_out = snakemake@output[['samples']]
comparison_entry = snakemake@params[['comparison']]
comparison_name = snakemake@params[['comparison_name']]

normalize_group_spec <- function(x) {
  # supports:
  # 1) named list: list(treatment="EtOH", time=c(1,2))
  # 2) list of named lists: list(list(treatment="EtOH"), list(time=c(1,2)))
  if (!is.list(x)) stop("Group spec must be a list.")

  # already named list of covariates
  if (!is.null(names(x)) && all(nzchar(names(x)))) return(x)

  # list-of-lists -> flatten to named covariate list
  out <- list()
  for (item in x) {
    if (!is.list(item) || is.null(names(item)) || length(item) != 1) {
      stop("Invalid group spec item; expected one named covariate per list item.")
    }
    nm <- names(item)[1]
    out[[nm]] <- item[[1]]
  }
  out
}

build_mask <- function(df, spec) {
  m <- rep(TRUE, nrow(df))
  for (covar in names(spec)) {
    if (!covar %in% colnames(df)) {
      stop(sprintf("Covariate '%s' not found in samplesheet.", covar))
    }
    vals <- spec[[covar]]
    # robust for factor columns + numeric/string values from YAML
    m <- m & (as.character(df[[covar]]) %in% as.character(vals))
  }
  m
}

# comparison_entry is a dictionary with entries:
# - type: twogroup
# - design: formula (optional)
# - groups:
#    - group1:
#      - factor1: value1, value2, ...
#    - group2:
#      - factor1: value1, value2, ...


# - [ ] multiple values -> inclusive or
# - [ ] no design given, fallback to simple additive of all (from samplesheet)                                                                                                                                                                              
# - [ ] inbalanced designs

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

# Dummy samplesheet
## Since designs can be inbalanced, we cannot just assume balance and use colMeans trick.
## Rather then some awkward balance check, let's just create a dummy samplesheet that is balanced for sure makes the trick robust for both balanced and inbalanced cases.
dummy_samplesheet <- unique(samplesheet)
dummy_samplesheet <- dummy_samplesheet[rep(seq_len(nrow(dummy_samplesheet)), each = 2), , drop = FALSE]

# Define contrast using the dummy samplesheet.
group_names <- setdiff(names(comparison_entry), c("type", "design"))
print(group_names)
group_specs <- lapply(group_names, function(g) normalize_group_spec(comparison_entry[[g]]))
names(group_specs) <- group_names
group_masks_dummy <- lapply(group_specs, function(spec) build_mask(dummy_samplesheet, spec))
group_masks_real <- lapply(group_specs, function(spec) build_mask(samplesheet, spec))
group_samples <- lapply(group_masks_real, function(m) rownames(samplesheet)[m])
relevantsamples <- unname(unique(unlist(group_samples, use.names = FALSE)))
print(paste("Relevant samples for", comparison_name, ":", paste(relevantsamples, collapse=", ")))
write_lines(relevantsamples,relevantsamples_out, sep='\n')

dummy_design <- model.matrix(formula, data = dummy_samplesheet)
c1 <- colMeans(dummy_design[group_masks_dummy[[1]], , drop = FALSE])
c2 <- colMeans(dummy_design[group_masks_dummy[[2]], , drop = FALSE])
contrast <- c2 - c1
print("Contrast set to: ", contrast)

# Run edgeR.
keep <- filterByExpr(
    countmat, design=design,min.count = 5, min.prop = 0.49
)
countmat <- countmat[keep,]
countmat <- countmat + pseudocount
countmat_disp <- estimateGLMCommonDisp(countmat, design, verbose=TRUE)
fit <- glmQLFit(countmat, design=design, dispersion = countmat_disp)
res <- glmQLFTest(fit, contrast=contrast)
results <- topTags(res, n = Inf, adjust.method="BH", sort.by='PValue', p.value=1)
res <- results$table
res$peak_id <- rownames(res)
res$group_assignment <- ifelse(res$logFC > 0, group_names[2], group_names[1])
write.table(
    res,
    file=edgeR_out,
    quote=FALSE,
    sep='\t',
    row.names=FALSE
)
