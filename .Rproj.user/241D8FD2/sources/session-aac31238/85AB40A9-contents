# ============================================================
# 04_ml_pipeline.R
# Module 3: Dimensionality reduction (PCA + UMAP),
#           K-means clustering, ElasticNet feature selection,
#           Random Forest response prediction
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(randomForest)
  library(cluster)
  library(pROC)
  library(purrr)
})


# ── 3.1  Dimensionality Reduction ────────────────────────────────────────────

dimensionality_reduction <- function(expr_filtered, demo) {
  cat("\n--- [3.1] Dimensionality Reduction ---\n")

  # PCA (centre + scale)
  pca_res   <- prcomp(expr_filtered, center = TRUE, scale. = TRUE)
  var_exp   <- cumsum(summary(pca_res)$importance[2, ])
  n90       <- which(var_exp >= 0.90)[1]
  cat(sprintf("  PCA: %d PCs explain 90%% of variance\n", n90))

  # UMAP via umap package if available, else t-SNE via Rtsne, else PCA 2D
  n_input  <- min(20, n90)
  pc_input <- pca_res$x[, seq_len(n_input)]

  use_umap  <- requireNamespace("umap",  quietly = TRUE)
  use_tsne  <- requireNamespace("Rtsne", quietly = TRUE)

  if (use_umap) {
    cat(sprintf("  Fitting UMAP on top %d PCs …\n", n_input))
    cfg               <- umap::umap.defaults
    cfg$n_neighbors   <- 15
    cfg$min_dist      <- 0.10
    cfg$random_state  <- 42
    emb <- umap::umap(pc_input, config = cfg)$layout
    dim_names <- c("UMAP1","UMAP2")
    method    <- "UMAP"
  } else if (use_tsne) {
    cat(sprintf("  Fitting t-SNE on top %d PCs …\n", n_input))
    set.seed(42)
    emb <- Rtsne::Rtsne(pc_input, dims = 2, perplexity = 30,
                         check_duplicates = FALSE, pca = FALSE)$Y
    dim_names <- c("tSNE1","tSNE2")
    method    <- "t-SNE"
  } else {
    cat("  [UMAP/Rtsne unavailable — using PC1/PC2 for 2D embedding]\n")
    emb       <- pca_res$x[, 1:2]
    dim_names <- c("PC1","PC2")
    method    <- "PCA"
  }

  embedding_df <- tibble::tibble(
    patient_id       = rownames(expr_filtered),
    !!dim_names[1]  := emb[, 1],
    !!dim_names[2]  := emb[, 2]
  ) |>
    dplyr::left_join(
      demo |> dplyr::select(patient_id, treatment, true_responder, treatment_label),
      by = "patient_id"
    )

  cat(sprintf("  %s embedding: %d patients\n", method, nrow(embedding_df)))

  list(pca_res = pca_res, embedding_df = embedding_df,
       pc_input = pc_input, method = method, dim_names = dim_names)
}


# ── 3.2  Patient Clustering ──────────────────────────────────────────────────

patient_clustering <- function(pca_res, demo, n_clusters = 3L) {
  cat(sprintf("\n--- [3.2] Patient Clustering (K-means, K=%d) ---\n", n_clusters))

  pc_mat <- pca_res$x[, 1:10]   # top 10 PCs

  # Elbow analysis (K = 2…6)
  set.seed(42)
  inertias <- vapply(2:6, function(k) {
    km <- kmeans(pc_mat, centers = k, nstart = 20, iter.max = 300)
    km$tot.withinss
  }, numeric(1))
  names(inertias) <- as.character(2:6)

  # Final K-means
  set.seed(42)
  km_fit <- kmeans(pc_mat, centers = n_clusters, nstart = 30, iter.max = 300)

  # Silhouette score
  sil_avg <- mean(cluster::silhouette(km_fit$cluster, dist(pc_mat))[, 3])
  cat(sprintf("  Silhouette score (K=%d): %.3f  [> 0.2 considered acceptable]\n",
              n_clusters, sil_avg))

  cluster_df <- tibble::tibble(
    patient_id = rownames(pca_res$x),
    cluster    = factor(km_fit$cluster)
  ) |>
    dplyr::left_join(
      demo |> dplyr::select(patient_id, true_responder, treatment),
      by = "patient_id"
    )

  cat("\n  Cluster composition (mean rates):\n")
  comp <- cluster_df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      N              = dplyr::n(),
      Responder_rate = mean(true_responder),
      Active_rate    = mean(treatment),
      .groups        = "drop"
    )
  print(comp)

  list(cluster_df = cluster_df, km_fit = km_fit, inertias = inertias)
}


# ── 3.3  Predictive Biomarker Modelling ──────────────────────────────────────

biomarker_response_prediction <- function(expr_filtered, proteomics, demo) {
  cat("\n--- [3.3] Predictive Biomarker Modelling ---\n")

  # Feature matrix: transcriptomics + baseline proteomics
  prot_bl <- proteomics |>
    dplyr::filter(week == 0) |>
    dplyr::select(-week, -treatment, -true_responder)

  protein_cols <- setdiff(names(prot_bl), "patient_id")

  common_ids <- intersect(rownames(expr_filtered), prot_bl$patient_id)

  X_trans <- expr_filtered[common_ids, , drop = FALSE]
  X_prot  <- prot_bl |>
    dplyr::filter(patient_id %in% common_ids) |>
    tibble::column_to_rownames("patient_id") |>
    as.matrix()
  X_prot  <- X_prot[common_ids, , drop = FALSE]

  X_all  <- cbind(X_trans, X_prot)
  y      <- demo$true_responder[match(common_ids, demo$patient_id)]

  cat(sprintf("  Feature matrix: %d patients \u00d7 %d features (%d genes + %d proteins)\n",
              nrow(X_all), ncol(X_all), ncol(X_trans), ncol(X_prot)))
  cat(sprintf("  Class balance: %d responders / %d non-responders\n",
              sum(y), sum(y == 0)))

  # Scale features
  X_scaled <- scale(X_all)

  # ── Stage 1: ElasticNet feature selection ─────────────────
  cat("\n  Stage 1: ElasticNet feature selection (cv.glmnet) …\n")
  set.seed(42)
  # Try multiple alpha values to find best mix
  alpha_vals <- c(0.3, 0.5, 0.7, 0.9)
  best_cvm   <- Inf
  best_enet  <- NULL
  best_alpha <- 0.5

  for (a in alpha_vals) {
    cv_fit <- cv.glmnet(X_scaled, y, alpha = a, family = "binomial",
                        nfolds = 5, standardize = FALSE, type.measure = "auc")
    if (min(cv_fit$cvm) < best_cvm) {
      best_cvm  <- min(cv_fit$cvm)
      best_enet <- cv_fit
      best_alpha <- a
    }
  }

  coef_vec    <- as.numeric(coef(best_enet, s = "lambda.min"))[-1]  # drop intercept
  selected    <- which(coef_vec != 0)
  sel_names   <- colnames(X_scaled)[selected]
  cat(sprintf("  ElasticNet selected: %d / %d features  (alpha=%.1f)\n",
              length(selected), ncol(X_scaled), best_alpha))
  if (length(sel_names))
    cat(sprintf("  Top features: %s\n",
                paste(head(sel_names, 8), collapse = ", ")))

  # Fall back to top 50 if ElasticNet selects too few
  if (length(selected) < 5) {
    cat("  [Expanding to top 50 by |coef|]\n")
    selected  <- order(abs(coef_vec), decreasing = TRUE)[1:50]
    sel_names <- colnames(X_scaled)[selected]
  }

  X_sel <- X_scaled[, selected, drop = FALSE]

  # ── Stage 2: Random Forest with 5-fold CV ─────────────────
  cat("\n  Stage 2: Random Forest (5-fold stratified CV) …\n")
  set.seed(42)
  folds <- caret::createFolds(y, k = 5, list = TRUE)   # use caret if available
  # Fallback: manual stratified folds
  make_folds <- function(y, k = 5) {
    idx_0 <- which(y == 0); idx_1 <- which(y == 1)
    f0 <- split(sample(idx_0), cut(seq_along(idx_0), k, labels = FALSE))
    f1 <- split(sample(idx_1), cut(seq_along(idx_1), k, labels = FALSE))
    lapply(seq_len(k), function(i) c(f0[[i]], f1[[i]]))
  }
  folds <- make_folds(y, k = 5)

  auc_scores <- vapply(seq_along(folds), function(i) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    rf_tmp <- randomForest::randomForest(
      x = X_sel[train_idx, ], y = factor(y[train_idx]),
      ntree = 300, mtry = max(1L, floor(sqrt(ncol(X_sel)))),
      classwt = c("0" = 1, "1" = sum(y == 0)/sum(y == 1))
    )
    probs <- predict(rf_tmp, X_sel[test_idx, ], type = "prob")[, "1"]
    # Simple AUC via rank correlation
    auc_val <- suppressMessages(pROC::auc(pROC::roc(y[test_idx], probs,
                                                     quiet = TRUE)))
    as.numeric(auc_val)
  }, numeric(1))

  cat(sprintf("  CV AUROC: %.3f \u00b1 %.3f\n", mean(auc_scores), sd(auc_scores)))

  # Full-data RF for importance + OOB predictions
  set.seed(42)
  rf_full <- randomForest::randomForest(
    x = X_sel, y = factor(y),
    ntree   = 300,
    mtry    = max(1L, floor(sqrt(ncol(X_sel)))),
    classwt = c("0" = 1, "1" = sum(y == 0)/sum(y == 1)),
    importance = TRUE
  )

  # OOB predicted probabilities for ROC
  oob_probs <- rf_full$votes[, "1"]
  roc_res   <- pROC::roc(y, oob_probs, quiet = TRUE)
  auc_oob   <- as.numeric(pROC::auc(roc_res))
  roc_df    <- tibble::tibble(
    fpr = 1 - rev(roc_res$specificities),
    tpr = rev(roc_res$sensitivities)
  )

  importance_df <- tibble::tibble(
    feature    = sel_names,
    importance = randomForest::importance(rf_full)[, "MeanDecreaseGini"]
  ) |>
    dplyr::arrange(dplyr::desc(importance))

  cat("\n  Top 10 predictive biomarkers:\n")
  print(head(importance_df, 10), n = 10, row.names = FALSE)

  list(
    rf_model      = rf_full,
    importance_df = importance_df,
    roc_df        = roc_df,
    auc_cv        = mean(auc_scores),
    auc_oob       = auc_oob,
    selected_feats = sel_names
  )
}


# ── 3.4  Visualisation ────────────────────────────────────────────────────────

plot_ml_results <- function(dr_res, cluster_res, ml_res,
                            output_dir = "outputs") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  embedding_df <- dr_res$embedding_df
  dim_names    <- dr_res$dim_names
  method       <- dr_res$method
  d1 <- dim_names[1]; d2 <- dim_names[2]

  cluster_df <- cluster_res$cluster_df
  emb <- dplyr::left_join(embedding_df,
                           cluster_df |> dplyr::select(patient_id, cluster),
                           by = "patient_id")

  resp_colors <- c("1" = "#d62728", "0" = "#aec7e8")
  arm_colors  <- c("Active" = "#d62728", "Placebo" = "#1f77b4")
  clust_pal   <- c("1" = "#e377c2", "2" = "#7f7f7f", "3" = "#bcbd22",
                   "4" = "#17becf", "5" = "#9467bd")

  # ── A: embedding by response ───────────────────────────────
  pA <- ggplot(emb, aes(.data[[d1]], .data[[d2]],
                         colour = factor(true_responder))) +
    geom_point(alpha = 0.80, size = 1.8) +
    scale_colour_manual(values = resp_colors,
                        labels = c("0" = "Non-Responder", "1" = "Responder")) +
    labs(title = sprintf("A.  %s: Responder Status", method),
         colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── B: embedding by treatment ──────────────────────────────
  pB <- ggplot(emb, aes(.data[[d1]], .data[[d2]], colour = treatment_label,
                         shape = treatment_label)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_colour_manual(values = arm_colors) +
    scale_shape_manual( values = c("Active" = 16, "Placebo" = 15)) +
    labs(title = sprintf("B.  %s: Treatment Arm", method),
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── C: Cluster overlay ────────────────────────────────────
  pC <- ggplot(emb, aes(.data[[d1]], .data[[d2]], colour = cluster)) +
    geom_point(alpha = 0.80, size = 1.8) +
    scale_colour_manual(values = clust_pal,
                        labels = paste("Cluster", levels(emb$cluster))) +
    labs(title = "C.  K-means Patient Clusters",
         colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── D: Elbow plot ─────────────────────────────────────────
  elbow_df <- tibble::tibble(
    K        = 2:6,
    inertia  = cluster_res$inertias
  )
  pD <- ggplot(elbow_df, aes(K, inertia)) +
    geom_line(colour = "#1f77b4", linewidth = 1.5) +
    geom_point(size = 3, colour = "#1f77b4") +
    geom_vline(xintercept = 3, colour = "#d62728",
               linetype = "dashed", linewidth = 1.2) +
    annotate("text", x = 3.15, y = max(elbow_df$inertia) * 0.95,
             label = "Selected K=3", colour = "#d62728", size = 3.2, hjust = 0) +
    scale_x_continuous(breaks = 2:6) +
    labs(title = "D.  K-Means Elbow Plot",
         x = "Number of Clusters (K)", y = "Within-cluster SS (inertia)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  # ── E: Feature importance ─────────────────────────────────
  top15 <- head(ml_res$importance_df, 15)
  known_prots <- c("IL6","TNF","IL10","IL17A","BAFF","APRIL","IgG_total",
                   "IgG1","IgG4","CRP","SAA","FcRn","IL21","CXCL13",
                   "IgA","IgM","IL4","IL13","IFNg","IL1B")
  top15 <- top15 |>
    dplyr::mutate(type = ifelse(feature %in% known_prots, "Protein", "Gene"))

  pE <- ggplot(top15,
               aes(x = importance,
                   y = reorder(feature, importance),
                   fill = type)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c("Protein" = "#d62728", "Gene" = "#1f77b4")) +
    labs(title = "E.  Top Predictive Biomarkers",
         subtitle = "Random Forest — Mean Decrease in Gini",
         x = "Feature Importance", y = NULL, fill = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── F: ROC curve ──────────────────────────────────────────
  auc_label <- sprintf("OOB AUROC = %.3f", ml_res$auc_oob)

  pF <- ggplot(ml_res$roc_df, aes(fpr, tpr)) +
    geom_area(fill = "#d62728", alpha = 0.12) +
    geom_line(colour = "#d62728", linewidth = 2) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40") +
    annotate("text", x = 0.6, y = 0.15, label = auc_label,
             size = 3.8, colour = "#d62728", fontface = "bold") +
    coord_equal() +
    labs(title = "F.  ROC Curve: Response Prediction",
         subtitle = "Out-of-bag prediction (Random Forest)",
         x = "False Positive Rate", y = "True Positive Rate") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  # ── Combine 6 panels ──────────────────────────────────────
  combined <- (pA + pB + pC) / (pD + pE + pF) +
    patchwork::plot_annotation(
      title = "Module 3 — ML Pipeline: UMAP \u00b7 Clustering \u00b7 Predictive Biomarkers",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )

  out_path <- file.path(output_dir, "03_ml_pipeline.png")
  ggsave(out_path, combined, width = 16, height = 11, dpi = 200)
  cat(sprintf("\n  \u2713 Saved: %s\n", out_path))
}
