# .libPaths(R.home("library"))
suppressMessages(library("readr"))
suppressMessages(library("edgeR"))

mat = snakemake@input[['mat']]
samplesheet = snakemake@params[['samplesheet']]
pseudocount = as.integer(snakemake@params[['pseudocount']])
of = snakemake@params[['outputfolder']]
edgeR_out = snakemake@output[['table']]
comparison_entry = snakemake@params[['comparison']]


# # Read in samplesheet.
# # rownames are samples, columns are factors.
# samplesheet = read.table(
#     samplesheet,
#     row.names = 1,
#     header=TRUE,
#     stringsAsFactors=TRUE)
# #print(samplesheet[,factors[1]])
# # read in count matrix.
# countmat <- read.csv(mat, sep='\t')
# rows <- paste(countmat[,1],countmat[,2], countmat[,3], sep='|')
# countmat <- countmat[-c(1:3)]
# rownames(countmat) <- rows
# # make sure the columns of count matrix are equal to the rows of the samplesheet.
# countmat <- countmat[, rownames(samplesheet)] 

# # Create model matrix design.
# design <- model.matrix(
#     as.formula(
#         paste(
#              '~',
#              paste(
#                 colnames(samplesheet),
#                 collapse=interactionstr
#              ) 
#         )
#     ),
#     data=samplesheet
# )
# # Take the relevant samples:
# relevantsamples <- c(
#     rownames(samplesheet)[eval(parse(text=gr1_subset))],
#     rownames(samplesheet)[eval(parse(text=gr2_subset))]
# )
# # Write out for later.
# write_lines(
#     relevantsamples,
#     paste0(of, '/samples.txt', sep='')
# )

# # Define a contrast.
# c1 <- colMeans(design[eval(parse(text=gr1_subset)),])
# c2 <- colMeans(design[eval(parse(text=gr2_subset)),])
# contrast <- c2 - c1 

# # Run edgeR.
# keep <- filterByExpr(
#     countmat, design=design,min.count = 5, min.prop = 0.49
# )
# countmat <- countmat[keep,]
# countmat <- countmat + pseudocount
# countmat_disp <- estimateGLMCommonDisp(countmat, design, verbose=TRUE)
# # 
# fit <- glmQLFit(countmat, design=design, dispersion = countmat_disp)
# res <- glmQLFTest(fit, contrast=contrast)
# results <- topTags(res, n = Inf, adjust.method="BH", sort.by='PValue', p.value=1)
# res <- results$table
# res$peak_id <- rownames(res)
# write.table(
#     res,
#     file=edgeRtable,
#     quote=FALSE,
#     sep='\t',
#     row.names=FALSE
# )
