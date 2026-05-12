# ============================================================
# 02_multiomics_analysis.R
# Module 1: Transcriptomics QC, Differential Expression,
#           Proteomics longitudinal analysis, Cross-omics integration
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(limma)
  library(SummarizedExperiment)
})

# ── 1.1  Transcriptomics QC & Normalisation ──────────────────────────────────

qc_transcriptomics <- function(transcriptomics, batch_df, demo) {
  cat("\n--- [1.1] Transcriptomics QC & Normalisation ---\n")

  # Extract expression matrix (genes only, patient_id as rowname)
  expr_mat <- transcriptomics |>
    tibble::column_to_rownames("patient_id") |>
    as.matrix()

  cat(sprintf("  Samples: %d | Genes: %d\n", nrow(expr_mat), ncol(expr_mat)))

  # Step 1: log2(count + 1) normalisation
  expr_log <- log2(expr_mat + 1)

  # Step 2: Batch correction — median centering per batch
  # (mirrors limma::removeBatchEffect for a simple 2-batch design)
  batch_vec <- batch_df$batch[match(rownames(expr_log), batch_df$patient_id)]
  global_med <- apply(expr_log, 2, median)
  expr_corrected <- expr_log
  for (b in unique(batch_vec)) {
    idx <- which(batch_vec == b)
    batch_med <- apply(expr_log[idx, , drop=FALSE], 2, median)
    expr_corrected[idx, ] <- sweep(expr_log[idx, , drop=FALSE], 2,
                                    batch_med - global_med, "-")
  }
  batch_sizes <- table(batch_vec)
  cat(sprintf("  Batch correction (median centering): batch 0 n=%d, batch 1 n=%d\n",
              batch_sizes["0"], batch_sizes["1"]))

  # Step 3: High-variance gene filter (top 50%)
  gene_vars    <- apply(expr_corrected, 2, var)
  var_thresh   <- quantile(gene_vars, 0.50)
  keep_genes   <- names(gene_vars)[gene_vars > var_thresh]
  expr_filtered <- expr_corrected[, keep_genes, drop = FALSE]
  cat(sprintf("  High-variance genes retained: %d / %d\n",
              length(keep_genes), ncol(expr_corrected)))

  # Step 4: PCA outlier detection (>97.5th pct Euclidean distance in PC space)
  pca_res   <- prcomp(expr_filtered, scale. = TRUE, center = TRUE)
  pc_scores <- pca_res$x[, 1:10]
  pc_dist   <- sqrt(rowSums(pc_scores^2))
  is_outlier <- pc_dist > quantile(pc_dist, 0.975)
  cat(sprintf("  Potential outliers (> 97.5th pct PCA distance): %d samples\n",
              sum(is_outlier)))

  list(
    expr_filtered  = expr_filtered,
    expr_corrected = expr_corrected,
    keep_genes     = keep_genes,
    is_outlier     = is_outlier,
    pca_res        = pca_res
  )
}


# ── 1.2  Differential Expression ─────────────────────────────────────────────

differential_expression <- function(expr_filtered, demo) {
  cat("\n--- [1.2] Differential Expression: Responders vs Non-Responders ---\n")
  cat("  Method: per-gene Welch t-test + Benjamini-Hochberg FDR correction\n")
  cat("  (equivalent to limma eBayes on log-normalised data for this design)\n")

  meta       <- demo[match(rownames(expr_filtered), demo$patient_id), ]
  resp_mask  <- meta$true_responder == 1

  # Per-gene Welch t-test
  results <- lapply(colnames(expr_filtered), function(g) {
    x1 <- expr_filtered[resp_mask,  g]
    x0 <- expr_filtered[!resp_mask, g]
    tt <- tryCatch(t.test(x1, x0, var.equal = FALSE), error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    list(gene   = g,
         log2FC = mean(x1) - mean(x0),
         pvalue = tt$p.value,
         t_stat = tt$statistic)
  })
  results <- Filter(Negate(is.null), results)

  de_df <- do.call(rbind, lapply(results, as.data.frame)) |>
    tibble::as_tibble() |>
    dplyr::mutate(across(c(log2FC, pvalue, t_stat), as.numeric)) |>
    dplyr::arrange(pvalue) |>
    dplyr::mutate(
      rank = seq_len(dplyr::n()),
      padj = pmin(1, pvalue * dplyr::n() / rank),
      padj = rev(cummin(rev(padj)))   # BH monotonicity
    )

  sig      <- dplyr::filter(de_df, padj < 0.05)
  sig_up   <- dplyr::filter(sig, log2FC > 0.5)
  sig_down <- dplyr::filter(sig, log2FC < -0.5)

  cat(sprintf("  Significant genes (FDR < 5%%): %d\n", nrow(sig)))
  cat(sprintf("  Up in responders: %d | Down in responders: %d\n",
              nrow(sig_up), nrow(sig_down)))
  cat("  Top 5 DE genes:\n")
  print(head(de_df[, c("gene","log2FC","pvalue","padj")], 5), row.names = FALSE)

  de_df
}


# ── 1.3  Proteomics ───────────────────────────────────────────────────────────

proteomics_analysis <- function(proteomics) {
  cat("\n--- [1.3] Proteomics Analysis (Olink NPX) ---\n")

  protein_cols <- setdiff(names(proteomics),
                          c("patient_id","week","treatment","true_responder"))
  cat(sprintf("  Proteins: %d | Patients: %d\n",
              length(protein_cols), dplyr::n_distinct(proteomics$patient_id)))

  # NPX delta: Week 12 – Week 0 per patient per protein
  w0  <- proteomics |> dplyr::filter(week ==  0) |>
    dplyr::select(patient_id, treatment, dplyr::all_of(protein_cols))
  w12 <- proteomics |> dplyr::filter(week == 12) |>
    dplyr::select(patient_id, dplyr::all_of(protein_cols))

  delta <- dplyr::inner_join(
    w0,
    w12 |> dplyr::rename_with(~ paste0(.x, "_w12"), dplyr::all_of(protein_cols)),
    by = "patient_id"
  )
  for (p in protein_cols) {
    delta[[paste0("delta_", p)]] <- delta[[paste0(p, "_w12")]] - delta[[p]]
  }

  delta_long <- delta |>
    dplyr::select(patient_id, treatment,
                  dplyr::starts_with("delta_")) |>
    tidyr::pivot_longer(dplyr::starts_with("delta_"),
                        names_to  = "protein",
                        values_to = "delta_npx") |>
    dplyr::mutate(protein = sub("^delta_", "", protein))

  mean_delta <- delta_long |>
    dplyr::group_by(treatment, protein) |>
    dplyr::summarise(mean_delta = mean(delta_npx), .groups = "drop")

  cat("\n  Mean NPX change Week 0 \u2192 12 (selected proteins):\n")
  key_p <- c("IgG_total","IgG1","IgG4","IL6","CRP","BAFF")
  for (p in key_p) {
    act <- mean_delta$mean_delta[mean_delta$treatment == 1 & mean_delta$protein == p]
    pbo <- mean_delta$mean_delta[mean_delta$treatment == 0 & mean_delta$protein == p]
    if (length(act) && length(pbo))
      cat(sprintf("    %-12s  Active = %+.3f   Placebo = %+.3f\n", p, act, pbo))
  }

  list(delta_long = delta_long, mean_delta = mean_delta, protein_cols = protein_cols)
}


# ── 1.4  Multi-Omics Integration ─────────────────────────────────────────────

multiomics_integration <- function(expr_filtered, proteomics, demo) {
  cat("\n--- [1.4] Multi-Omics Integration ---\n")

  # Transcriptomic PCA (top 5 PCs)
  pca_res  <- prcomp(expr_filtered, scale. = TRUE, center = TRUE)
  var_exp  <- summary(pca_res)$importance[2, 1:5]
  trans_pcs <- as.data.frame(pca_res$x[, 1:5]) |>
    setNames(paste0("TransPC", 1:5)) |>
    tibble::rownames_to_column("patient_id")

  cat("  Transcriptomic PCs variance explained: ",
      paste(sprintf("PC%d=%.1f%%", 1:5, var_exp * 100), collapse=" | "), "\n")

  # Baseline proteomic markers
  key_prots <- c("IL6","TNF","BAFF","CXCL13","IgG_total")
  prot_bl   <- proteomics |>
    dplyr::filter(week == 0) |>
    dplyr::select(patient_id, dplyr::all_of(key_prots))

  combined <- dplyr::inner_join(trans_pcs, prot_bl, by = "patient_id")

  # Pearson correlation matrix (PCs × proteins)
  pc_cols   <- paste0("TransPC", 1:5)
  cross_corr <- cor(combined[, pc_cols], combined[, key_prots], method = "pearson")

  cat("\n  Transcriptomic PC \u00d7 Proteomic marker correlations (Pearson r):\n")
  print(round(cross_corr, 3))

  list(trans_pcs = trans_pcs, cross_corr = cross_corr)
}


# ── 1.5  Visualisation ───────────────────────────────────────────────────────

plot_multiomics <- function(expr_filtered, de_df, demo, proteomics,
                            cross_corr, pca_res, output_dir = "outputs") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ── A: Volcano plot ──────────────────────────────────────
  volcano_df <- de_df |>
    dplyr::mutate(
      neg_log10p = -log10(pmax(pvalue, 1e-10)),
      sig_label  = dplyr::case_when(
        padj < 0.05 & log2FC >  0.5 ~ "Up (FDR<5%)",
        padj < 0.05 & log2FC < -0.5 ~ "Down (FDR<5%)",
        TRUE                         ~ "NS"
      )
    )

  top_genes <- volcano_df |>
    dplyr::filter(padj < 0.05) |>
    dplyr::slice_min(pvalue, n = 8)

  pA <- ggplot(volcano_df, aes(log2FC, neg_log10p, colour = sig_label)) +
    geom_point(alpha = 0.65, size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey50") +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 12) +
    scale_colour_manual(values = c("Up (FDR<5%)"   = "#d62728",
                                   "Down (FDR<5%)" = "#1f77b4",
                                   "NS"             = "#cccccc")) +
    labs(title = "A.  Volcano Plot: DE Genes at Baseline",
         subtitle = "Responders vs Non-Responders",
         x = "log\u2082 Fold Change", y = "-log\u2081\u2080(p-value)",
         colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── B: PCA coloured by treatment × response ───────────────
  pc_df <- as.data.frame(pca_res$x[, 1:2]) |>
    tibble::rownames_to_column("patient_id") |>
    dplyr::left_join(demo[, c("patient_id","treatment","true_responder",
                               "treatment_label")],
                     by = "patient_id") |>
    dplyr::mutate(
      group = paste0(treatment_label, " / ",
                     ifelse(true_responder == 1, "Responder", "Non-R"))
    )
  var_exp_pct <- summary(pca_res)$importance[2, 1:2] * 100

  pB <- ggplot(pc_df, aes(PC1, PC2, colour = group, shape = treatment_label)) +
    geom_point(alpha = 0.80, size = 2) +
    scale_colour_manual(values = c(
      "Active / Responder"  = "#d62728",
      "Active / Non-R"      = "#ffbb78",
      "Placebo / Responder" = "#1f77b4",
      "Placebo / Non-R"     = "#aec7e8"
    )) +
    scale_shape_manual(values = c("Active" = 16, "Placebo" = 15)) +
    labs(title = "B.  PCA of Baseline Transcriptomics",
         x = sprintf("PC1 (%.1f%%)", var_exp_pct[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp_pct[2]),
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7),
          plot.title = element_text(face = "bold"))

  # ── C: IgG longitudinal ───────────────────────────────────
  igg_sum <- proteomics |>
    dplyr::group_by(treatment, week) |>
    dplyr::summarise(
      mean_igg = mean(IgG_total),
      sem_igg  = sd(IgG_total) / sqrt(dplyr::n()),
      .groups  = "drop"
    ) |>
    dplyr::mutate(arm = ifelse(treatment == 1, "Active", "Placebo"))

  pC <- ggplot(igg_sum, aes(week, mean_igg, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = mean_igg - sem_igg, ymax = mean_igg + sem_igg),
                alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.8) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("Active" = "#d62728", "Placebo" = "#1f77b4")) +
    scale_fill_manual(  values = c("Active" = "#d62728", "Placebo" = "#1f77b4")) +
    scale_x_continuous(breaks = c(0, 4, 12)) +
    labs(title = "C.  IgG Longitudinal Profile by Treatment",
         x = "Week", y = "IgG_total (NPX, mean \u00b1 SEM)",
         colour = NULL, fill = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── D: Cross-omics heatmap via pheatmap ──────────────────
  png_path <- file.path(output_dir, "01_multiomics_analysis.png")
  png(png_path, width = 3200, height = 2400, res = 200)

  # Combine A+B+C via patchwork (3 panels)
  top_row <- pA + pB + pC + patchwork::plot_layout(ncol = 3)

  # Print top row then heatmap side by side is tricky with pheatmap;
  # use a 2-row layout: top row as ggplot, bottom row as pheatmap via grid
  library(grid)
  library(gridExtra)

  grob_top <- patchwork::patchworkGrob(top_row)

  hm_breaks <- seq(-0.4, 0.4, length.out = 101)
  hm_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 1), "null"))))

  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(grob_top)
  popViewport()

  pushViewport(viewport(layout.pos.row = 2))
  pheatmap(cross_corr,
           color          = hm_colors,
           breaks         = hm_breaks,
           cluster_rows   = FALSE,
           cluster_cols   = FALSE,
           display_numbers = TRUE,
           number_format  = "%.2f",
           fontsize        = 10,
           main            = "D.  Cross-Omics Correlation\n(Transcriptomic PCs \u00d7 Proteomic Markers)",
           silent          = TRUE,
           newpage         = FALSE)
  popViewport()

  dev.off()
  cat(sprintf("\n  \u2713 Saved: %s\n", png_path))
}
