
debug_flag <- FALSE
# check if DEBUG flag is set
if (!is.null(snakemake@config$DEBUG)) {
  message("debug flag is set")
  # if set, then check if True
  if (snakemake@config$DEBUG) {
    debug_flag <- TRUE
    message("In debug mode: saving R objects to inspect later")
    path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
    dir.create(path_debug, showWarnings = FALSE)
    save(snakemake, file = file.path(path_debug, "doublet_scores.rdata"))
  }
}

merged = data.frame(matrix(ncol=2, nrow=0))

for (i in 1:length(snakemake@params$samples)){
	temp = read.csv(snakemake@input[[i]], stringsAsFactors= FALSE)
	temp$cell = paste0(snakemake@params$samples[[i]], '_',temp$cell)
	merged = rbind(merged, temp)
}

write.csv(merged, file = snakemake@output[[1]], quote=FALSE, row.names=FALSE)


if (debug_flag) {
  save.image(file = file.path(path_debug, "doublet_scores_workspace.rdata"))
}
