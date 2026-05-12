# ============================================================
# 01_simulate_data.R
# Simulate a Phase 2 autoimmune trial with an FcRn inhibitor
# (VYVGART-like) in a severe autoimmune disease (e.g. gMG).
#
# Data modalities:
#   - Demographics & treatment assignment
#   - Longitudinal IgG + disease score (weeks 0,4,8,12,24)
#   - Baseline transcriptomics (500 genes, RNA-seq-like)
#   - Longitudinal Olink proteomics (50 proteins, weeks 0,4,12)
#   - Time-to-response for survival analysis
# ============================================================

simulate_trial_data <- function(n_patients  = 120,
                                n_genes     = 500,
                                n_proteins  = 50,
                                timepoints  = c(0, 4, 8, 12, 24),
                                seed        = 42) {
  set.seed(seed)

  patient_ids <- sprintf("PT%03d", seq_len(n_patients))

  # ── Treatment (1:1 randomisation) ────────────────────────
  treatment <- rbinom(n_patients, 1, 0.5)

  # ── Demographics ─────────────────────────────────────────
  age              <- pmax(18, pmin(75, rnorm(n_patients, 45, 12)))
  sex              <- rbinom(n_patients, 1, 0.65)          # 65% female
  disease_duration <- pmax(0.5, rexp(n_patients, rate = 0.2))
  baseline_igg     <- pmax(5,   pmin(25, rnorm(n_patients, 12, 3)))

  # ── Latent biology → responder status ────────────────────
  # Active ~65% responders; Placebo ~20%
  latent    <- rnorm(n_patients)
  threshold <- ifelse(treatment == 1, -0.5, 0.8)
  responder <- as.integer(latent > threshold)

  demo <- tibble::tibble(
    patient_id        = patient_ids,
    treatment         = treatment,
    treatment_label   = ifelse(treatment == 1, "Active", "Placebo"),
    age               = round(age, 1),
    sex               = sex,
    sex_label         = ifelse(sex == 1, "F", "M"),
    disease_duration  = round(disease_duration, 2),
    baseline_igg      = round(baseline_igg, 2),
    true_responder    = responder,
    latent_biology    = round(latent, 4)
  )

  # ── Longitudinal (IgG + disease score) ───────────────────
  long_list <- vector("list", n_patients * length(timepoints))
  idx <- 1L
  for (i in seq_len(n_patients)) {
    for (wk in timepoints) {
      trt  <- treatment[i]
      resp <- responder[i]
      bigg <- baseline_igg[i]

      # IgG PD: Emax-style reduction
      if (trt == 1) {
        max_red    <- ifelse(resp == 1, 0.65, 0.30)
        igg_factor <- 1 - max_red * min(wk / 8, 1) * rnorm(1, 1, 0.08)
      } else {
        igg_factor <- 1 + rnorm(1, 0, 0.04)
      }
      igg <- pmax(1, pmin(25, bigg * igg_factor + rnorm(1, 0, 0.4)))

      # Disease score 0-10 (lower = better)
      base_score <- rnorm(1, 6.5, 1.5)
      improve    <- dplyr::case_when(
        trt == 1 & resp == 1 ~ 3.0 * min(wk / 12, 1),
        trt == 1             ~ 0.8 * min(wk / 12, 1),
        TRUE                 ~ 0.2 * min(wk / 12, 1)
      )
      score <- pmax(0, pmin(10, base_score - improve + rnorm(1, 0, 0.45)))

      long_list[[idx]] <- tibble::tibble(
        patient_id     = patient_ids[i],
        week           = wk,
        treatment      = trt,
        true_responder = resp,
        igg_gL         = round(igg, 2),
        disease_score  = round(score, 2),
        pct_igg_change = round((igg - bigg) / bigg * 100, 1)
      )
      idx <- idx + 1L
    }
  }
  longitudinal <- dplyr::bind_rows(long_list)

  # ── Transcriptomics (log-normal counts) ──────────────────
  gene_names <- sprintf("GENE%04d", seq_len(n_genes))
  expr_mat   <- matrix(rlnorm(n_patients * n_genes, 3, 1.5),
                       nrow = n_patients, ncol = n_genes,
                       dimnames = list(patient_ids, gene_names))

  # 50 DE genes up in responders at baseline
  de_idx <- sample(n_genes, 50, replace = FALSE)
  for (idx in de_idx) {
    resp_mask          <- responder == 1
    expr_mat[resp_mask, idx] <- expr_mat[resp_mask, idx] *
      runif(sum(resp_mask), 1.5, 3.0)
  }

  # Batch effect (2 sequencing runs)
  batch        <- rbinom(n_patients, 1, 0.5)
  batch_effect <- rnorm(n_genes, 0, 0.3)
  for (i in which(batch == 1)) {
    expr_mat[i, ] <- expr_mat[i, ] * exp(batch_effect)
  }

  transcriptomics <- as.data.frame(expr_mat) |>
    tibble::rownames_to_column("patient_id")

  batch_df <- tibble::tibble(patient_id = patient_ids, batch = batch)

  # ── Proteomics (Olink NPX, 3 timepoints) ─────────────────
  named_prots <- c("IL6","TNF","IL10","IL17A","BAFF","APRIL",
                   "IgG_total","IgG1","IgG4","CRP","SAA","FcRn",
                   "IL21","CXCL13","IgA","IgM","IL4","IL13","IFNg","IL1B")
  anon_prots  <- sprintf("PROT%03d", seq_len(n_proteins - length(named_prots)))
  prot_names  <- c(named_prots, anon_prots)

  prot_list <- vector("list", n_patients * 3)
  idx <- 1L
  for (i in seq_len(n_patients)) {
    base_prot <- rnorm(n_proteins, 5, 2)
    if (responder[i] == 1)
      base_prot[c(1,2,3,4,5,6,10,11)] <-
        base_prot[c(1,2,3,4,5,6,10,11)] + runif(8, 0.5, 1.5)

    for (wk in c(0, 4, 12)) {
      prot <- base_prot
      if (treatment[i] == 1) {
        sc      <- min(wk / 8, 1)
        sc12    <- min(wk / 12, 1)
        emax_r  <- ifelse(responder[i] == 1, 2.0, 1.0)
        prot[7] <- prot[7] - emax_r * sc           # IgG_total
        prot[8] <- prot[8] - (emax_r + 0.5) * sc   # IgG1
        prot[9] <- prot[9] - (emax_r - 0.5) * sc   # IgG4
        if (responder[i] == 1)
          prot[c(1,2,10,11)] <- prot[c(1,2,10,11)] - 0.8 * sc12
      }
      prot <- prot + rnorm(n_proteins, 0, 0.3)

      row_df <- tibble::tibble(
        patient_id     = patient_ids[i],
        week           = wk,
        treatment      = treatment[i],
        true_responder = responder[i]
      )
      vals <- as.list(setNames(round(prot, 3), prot_names))
      prot_list[[idx]] <- dplyr::bind_cols(row_df, tibble::as_tibble(vals))
      idx <- idx + 1L
    }
  }
  proteomics <- dplyr::bind_rows(prot_list)

  # ── Survival (time-to-response) ───────────────────────────
  surv_list <- vector("list", n_patients)
  for (i in seq_len(n_patients)) {
    if (responder[i] == 1) {
      tte   <- if (treatment[i] == 1)
                 rweibull(1, shape = 2, scale = 9) + 2
               else
                 rweibull(1, shape = 1.5, scale = 15) + 5
      event <- 1L
    } else {
      tte   <- runif(1, 20, 24)
      event <- 0L
    }
    # ~10% censoring
    if (runif(1) < 0.10) {
      tte   <- runif(1, 4, tte)
      event <- 0L
    }
    surv_list[[i]] <- tibble::tibble(
      patient_id             = patient_ids[i],
      treatment              = treatment[i],
      true_responder         = responder[i],
      time_to_response_weeks = round(pmin(tte, 24), 1),
      event_observed         = event,
      baseline_igg           = baseline_igg[i],
      latent_biology         = latent[i]
    )
  }
  survival <- dplyr::bind_rows(surv_list)

  list(
    demographics    = demo,
    longitudinal    = longitudinal,
    transcriptomics = transcriptomics,
    batch           = batch_df,
    proteomics      = proteomics,
    survival        = survival,
    de_gene_indices = de_idx,
    gene_names      = gene_names,
    protein_names   = prot_names
  )
}
