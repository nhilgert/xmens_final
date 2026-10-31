# Reproduce Monocle 2's plot_pseudotime_heatmap() on a cell subset (default: the 3 dpi panel).
# (plot_pseudotime_heatmap was removed from Monocle 3; papers that cite Monocle 3 and show this
#  heatmap either call the Monocle 2 function or re-implement it.)
# Ported from cole-trapnell-lab/monocle-release R/plotting.R, R/expr_models.R, R/utils.R:
#   size factors   : estimateSizeFactorsForMatrix(method = "mean-geometric-mean-total")
#   dispersion     : disp_calc_helper_NB + parametricDispersionFit (+ Cook's-distance outlier refit)
#   smoother       : fit_model_helper -> VGAM::vglm(round(x/SF) ~ sm.ns(Pseudotime, df=3),
#                    negbinomial.size(size = 1/disp_func(mean(round(x)))), epsilon = 1e-1)
#   curve          : predict(type="response") on 100 evenly spaced Pseudotime points
#   transform      : log10(m + 1); row z-score; clip to [-3, 3]
#   drawing        : pheatmap, breaks seq(-3.1, 3.1, 0.1), colorRamps::blue2green2red(61),
#                    rows clustered (ward.D2 on (1 - cor)/2), cutree 6
# usage: Rscript monocle_heatmap_repro.R PERCELL_TSV_GZ GENE_ORDER_TSV OUTDIR [SUBSET] [DISP_TABLE_TSV_GZ]
#   SUBSET = a timepoint label (default 3_dpi), 'all' for every cell on the global axis, or 'induced' for every
#            timepoint except uninduced (pseudotime re-ranked within those cells; the pseudotime itself is not refit)
#   DISP_TABLE = optional per-gene (gene, mu, disp) moments for every gene detected in >1 cell of the same object
#            (e.g. disp_table_mouse.tsv.gz); the dispersion trend is then fit genome-wide, as Monocle's
#            estimateDispersions does, instead of over the panel genes only
suppressPackageStartupMessages({library(VGAM); library(pheatmap); library(colorRamps); library(splines)})
args <- commandArgs(trailingOnly = TRUE)
percell <- args[1]; order_tsv <- args[2]; outdir <- args[3]
subset_lab <- if (length(args) >= 4) args[4] else "3_dpi"
disp_tsv <- if (length(args) >= 5) args[5] else NA
dir.create(outdir, showWarnings = FALSE)

pc <- read.delim(gzfile(percell), check.names = FALSE)
genes <- read.delim(order_tsv)$gene
cnt <- as.matrix(pc[, paste0("cnt_", genes)]); colnames(cnt) <- genes

# ---- size factors over the whole object (Monocle stores them on the full cds) -----------------
sf_all <- pc$total_counts / exp(mean(log(pc$total_counts)))

# ---- dispersion trend (estimateDispersions on the full cds; here only the 28 panel genes exist) -
disp_table <- local({
  x <- t(round(t(cnt)) )                       # cells x genes, already integer
  xn <- sweep(x, 1, sf_all, "/")
  mu <- colMeans(xn); v <- colMeans(sweep(xn, 2, mu)^2)
  xim <- mean(1 / sf_all)
  d <- (v - xim * mu) / mu^2; d[d < 0] <- 0
  data.frame(mu = mu, disp = d, row.names = genes)
})
parametricDispersionFit <- function(disp_table, initial_coefs = c(1e-6, 1)) {
  coefs <- initial_coefs; iter <- 0
  while (TRUE) {
    residuals <- disp_table$disp / (coefs[1] + coefs[2] / disp_table$mu)
    good <- disp_table[which((residuals > initial_coefs[1]) & (residuals < 10000)), ]
    fit <- suppressWarnings(glm(disp ~ I(1/mu), data = good, family = Gamma(link = "identity"), start = coefs))
    oldcoefs <- coefs; coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]) coefs[1] <- initial_coefs[1]
    if (coefs[2] < 0) stop("Parametric dispersion fit failed")
    if (sum(log(coefs / oldcoefs)^2) < initial_coefs[1]) break
    iter <- iter + 1; if (iter > 10) { warning("Dispersion fit did not converge."); break }
  }
  list(fit, coefs)
}
if (!is.na(disp_tsv)) {                       # genome-wide moments; clamp as Monocle's disp_calc_helper_NB does
  disp_table <- read.delim(gzfile(disp_tsv))
  rownames(disp_table) <- make.unique(disp_table$gene)
  disp_table <- disp_table[is.finite(disp_table$mu) & disp_table$mu > 0, c("mu", "disp")]
  disp_table$disp[!is.finite(disp_table$disp) | disp_table$disp < 0] <- 0
}
res <- parametricDispersionFit(disp_table)
CD <- cooks.distance(res[[1]]); outl <- union(names(CD[CD > 4 / nrow(disp_table)]), setdiff(rownames(disp_table), names(CD)))
res <- parametricDispersionFit(disp_table[!rownames(disp_table) %in% outl, ])
coefs <- res[[2]]
disp_func <- function(q) coefs[1] + coefs[2] / q
cat(sprintf("dispersion trend over %d genes: asymptDisp=%.4g extraPois=%.4g (outliers removed: %d)\n", nrow(disp_table), coefs[1], coefs[2], length(outl)))

# ---- cell subset ----------------------------------------------------------------------------
keep <- if (subset_lab == "all") rep(TRUE, nrow(pc)) else if (subset_lab == "induced") pc$timepoint != "uninduced" else pc$timepoint == subset_lab
d3 <- pc[keep, ]; c3 <- cnt[keep, , drop = FALSE]; sf3 <- sf_all[keep]
pt_inputs <- if (subset_lab == "all") list(aligned = pc$pt_rank_norm) else if (subset_lab == "induced")
  list(aligned = (rank(d3$pt_rank_norm, ties.method = "first") - 1) / (sum(keep) - 1)) else list(
  aligned = (rank(d3$pt_rank_norm, ties.method = "first") - 1) / (sum(keep) - 1),  # within-3 dpi rank, 0..1
  global  = d3$pt_rank_norm                                                          # global rank-norm value as given
)

smooth_gene <- function(x, Pseudotime, newdata, family = c("trend", "free")) {
  family <- match.arg(family)
  f_expression <- round(x / sf3)
  fam <- if (family == "trend") negbinomial.size(size = 1 / disp_func(mean(round(x)))) else negbinomial()
  fit <- tryCatch(suppressWarnings(vglm(f_expression ~ sm.ns(Pseudotime, df = 3), family = fam, epsilon = 1e-1)),
                  error = function(e) suppressWarnings(vglm(f_expression ~ sm.ns(Pseudotime, df = 3),
                                                            family = negbinomial(), epsilon = 1e-1, checkwz = TRUE)))
  as.numeric(predict(fit, newdata = newdata, type = "response"))
}
to_heat <- function(m, scale_max = 3, scale_min = -3) {
  m <- m[!apply(m, 1, sum) == 0, , drop = FALSE]
  m <- log10(m + 1)
  m <- m[!apply(m, 1, sd) == 0, , drop = FALSE]
  m <- t(scale(t(m), center = TRUE))
  m[is.nan(m)] <- 0; m[m > scale_max] <- scale_max; m[m < scale_min] <- scale_min
  m
}

for (nm in names(pt_inputs)) {
  Pseudotime <- pt_inputs[[nm]]
  newdata <- data.frame(Pseudotime = seq(min(Pseudotime), max(Pseudotime), length.out = 100))
  for (fam in c("trend", "free")) {
    m <- t(sapply(genes, function(g) smooth_gene(c3[, g], Pseudotime, newdata, fam)))
    write.table(data.frame(gene = genes, m, check.names = FALSE), file.path(outdir, sprintf("monocle_curve_%s_%s.tsv", nm, fam)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data.frame(gene = rownames(to_heat(m)), to_heat(m), check.names = FALSE),
                file.path(outdir, sprintf("monocle_heat_%s_%s.tsv", nm, fam)), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write.table(newdata, file.path(outdir, sprintf("monocle_grid_%s.tsv", nm)), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- palette and a faithful pheatmap render (aligned axis, dispersion trend) -------------------
bks <- seq(-3.1, 3.1, by = 0.1); hmcols <- blue2green2red(length(bks) - 1)
writeLines(hmcols, file.path(outdir, "blue2green2red_61.txt"))
hm <- as.matrix(read.delim(file.path(outdir, "monocle_heat_aligned_trend.tsv"), row.names = 1, check.names = FALSE))
colnames(hm) <- 1:ncol(hm)
row_dist <- as.dist((1 - cor(t(hm))) / 2); row_dist[is.na(row_dist)] <- 1
ph <- pheatmap(hm, useRaster = TRUE, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE,
               clustering_distance_rows = row_dist, clustering_method = "ward.D2", cutree_rows = 6, silent = TRUE,
               filename = NA, breaks = bks, border_color = NA, color = hmcols)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 6)))
for (clust in c(TRUE, FALSE)) {
  ph_res <- pheatmap(hm, useRaster = TRUE, cluster_cols = FALSE, cluster_rows = clust, show_rownames = TRUE, show_colnames = FALSE,
                     clustering_distance_rows = row_dist, clustering_method = "ward.D2", cutree_rows = if (clust) 6 else NA,
                     annotation_row = if (clust) annotation_row else NA, treeheight_row = 20, breaks = bks, fontsize = 6,
                     color = hmcols, border_color = NA, silent = TRUE, filename = NA)
  f <- file.path(outdir, sprintf("monocle_pheatmap_%s_%s.png", gsub("_", "", subset_lab), if (clust) "clustered" else "published_order"))
  png(f, width = 4, height = 5, units = "in", res = 300); grid::grid.newpage(); grid::grid.draw(ph_res$gtable); dev.off()
  if (clust) writeLines(rownames(hm)[ph$tree_row$order], file.path(outdir, "monocle_clustered_row_order.txt"))
}
cat("done\n")
