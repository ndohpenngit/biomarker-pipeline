# ============================================================
# 03_longitudinal_survival.R
# Module 2: Linear Mixed-Effects Model, Emax PD model,
#           Kaplan-Meier, Cox Proportional Hazards
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(survival)
  library(broom)
})


# ── 2.1  Linear Mixed-Effects Model ──────────────────────────────────────────

fit_mixed_effects_model <- function(longitudinal) {
  cat("\n--- [2.1] Linear Mixed-Effects Model ---\n")
  cat("  Spec: disease_score ~ week * treatment + (1 + week | patient_id)\n")

  # lme4::lmer with lmerTest for Satterthwaite df and p-values
  lme_fit <- lme4::lmer(
    disease_score ~ week * treatment + (1 + week | patient_id),
    data    = longitudinal,
    REML    = TRUE,
    control = lme4::lmerControl(optimizer = "bobyqa",
                                 optCtrl  = list(maxfun = 2e5))
  )

  cat("\n  Fixed-Effect Estimates:\n")
  # Use lmerTest::summary if available, otherwise lme4 summary
  fe_summary <- tryCatch(
    coef(lmerTest::as_lmerModLmerTest(lme_fit))[["time_fixed"]],
    error = function(e) summary(lme_fit)$coefficients
  )
  fe_summary <- summary(lme_fit)$coefficients
  print(round(fe_summary, 4))

  # Drug effect at week 24
  fe     <- lme4::fixef(lme_fit)
  eff_24 <- fe["week:treatment"] * 24
  cat(sprintf("\n  Estimated cumulative drug effect at Week 24: %+.2f score units\n", eff_24))
  cat("  (Negative = Active arm improves more than Placebo)\n")

  lme_fit
}


# ── 2.2  Emax Pharmacodynamic Model ──────────────────────────────────────────

fit_emax_pd_model <- function(longitudinal) {
  cat("\n--- [2.2] Emax Pharmacodynamic Model (IgG reduction, Active arm) ---\n")

  emax_fn <- function(week, emax, ec50) -emax * week / (ec50 + week)

  active_df <- longitudinal |> dplyr::filter(treatment == 1)

  params_list <- list()
  for (resp in c(1L, 0L)) {
    label  <- ifelse(resp == 1, "Responders", "Non-Responders")
    sub_df <- active_df |>
      dplyr::filter(true_responder == resp) |>
      dplyr::group_by(week) |>
      dplyr::summarise(mean_pct = mean(pct_igg_change), .groups = "drop")

    tryCatch({
      fit <- nls(mean_pct ~ emax_fn(week, emax, ec50),
                 data  = sub_df,
                 start = list(emax = 50, ec50 = 4),
                 lower = c(0, 0.5), upper = c(100, 24),
                 algorithm = "port",
                 control   = nls.control(maxiter = 200))
      pars <- coef(fit)
      cat(sprintf("  %-20s  Emax = %.1f%%   EC50 = %.1f weeks\n",
                  label, pars["emax"], pars["ec50"]))
      params_list[[label]] <- pars
    }, error = function(e) {
      obs_max <- abs(min(sub_df$mean_pct))
      cat(sprintf("  %-20s  Observed max reduction \u2248 %.1f%%\n", label, obs_max))
      params_list[[label]] <<- c(emax = obs_max, ec50 = 6)
    })
  }

  list(active_df = active_df, params = params_list)
}


# ── 2.3  Kaplan-Meier Analysis ────────────────────────────────────────────────

kaplan_meier_analysis <- function(survival_df) {
  cat("\n--- [2.3] Kaplan-Meier Time-to-Response ---\n")

  # Fit KM curves per treatment arm
  km_fit <- survfit(
    Surv(time_to_response_weeks, event_observed) ~ treatment_label,
    data = survival_df |>
      dplyr::mutate(treatment_label = ifelse(treatment == 1, "Active", "Placebo"))
  )

  km_summary <- summary(km_fit)$table
  print(km_summary[, c("records","events","median","0.95LCL","0.95UCL")])

  # Log-rank test
  lr_test <- survdiff(
    Surv(time_to_response_weeks, event_observed) ~ treatment,
    data = survival_df
  )
  p_val <- 1 - pchisq(lr_test$chisq, df = 1)
  cat(sprintf("\n  Log-rank \u03c7\u00b2 = %.3f  p = %.5f\n", lr_test$chisq, p_val))

  list(km_fit = km_fit, lr_chisq = lr_test$chisq, lr_p = p_val)
}


# ── 2.4  Cox Proportional Hazards ────────────────────────────────────────────

cox_ph_analysis <- function(survival_df) {
  cat("\n--- [2.4] Cox Proportional Hazards Model ---\n")

  surv_df2 <- survival_df |>
    dplyr::mutate(
      baseline_igg_z    = scale(baseline_igg)[, 1],
      latent_biology_z  = scale(latent_biology)[, 1]
    )

  cox_fit <- coxph(
    Surv(time_to_response_weeks, event_observed) ~
      treatment + baseline_igg_z + latent_biology_z,
    data = surv_df2, ties = "efron"
  )

  cat("\n  Cox Model — Hazard Ratios (HR > 1 = faster time-to-response):\n")
  hr_df <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) |>
    dplyr::select(term, estimate, conf.low, conf.high, p.value) |>
    dplyr::rename(HR = estimate, `HR_lower` = conf.low,
                  `HR_upper` = conf.high, `p_value` = p.value)
  print(as.data.frame(round(hr_df[, -1], 3)), row.names = FALSE)

  trt_hr <- exp(coef(cox_fit)["treatment"])
  cat(sprintf("\n  treatment HR = %.2f\u00d7 higher instantaneous response rate in Active arm\n",
              trt_hr))

  list(cox_fit = cox_fit, hr_df = hr_df, trt_hr = trt_hr)
}


# ── 2.5  Visualisation ────────────────────────────────────────────────────────

plot_longitudinal_survival <- function(longitudinal, survival_df,
                                       lme_fit, emax_res, km_res, cox_res,
                                       output_dir = "outputs") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  arm_colors <- c("Active" = "#d62728", "Placebo" = "#1f77b4")

  # ── A: Disease score spaghetti + mean ─────────────────────
  long_plot <- longitudinal |>
    dplyr::mutate(arm = ifelse(treatment == 1, "Active", "Placebo"))

  mean_score <- long_plot |>
    dplyr::group_by(arm, week) |>
    dplyr::summarise(
      mean_s = mean(disease_score),
      sem_s  = sd(disease_score) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # Sample 20 patients per arm for spaghetti
  samp_ids <- long_plot |>
    dplyr::distinct(patient_id, arm) |>
    dplyr::group_by(arm) |>
    dplyr::slice_sample(n = 18) |>
    dplyr::pull(patient_id)

  samp_data <- long_plot |> dplyr::filter(patient_id %in% samp_ids)

  pA <- ggplot(mean_score, aes(week, mean_s, colour = arm, fill = arm)) +
    geom_line(data = samp_data,
              mapping = aes(x = week, y = disease_score,
                            colour = arm, group = patient_id),
              alpha = 0.10, linewidth = 0.5) +
    geom_ribbon(aes(ymin = mean_s - sem_s, ymax = mean_s + sem_s),
                alpha = 0.20, colour = NA) +
    geom_line(linewidth = 2.2) +
    geom_point(size = 3) +
    scale_colour_manual(values = arm_colors) +
    scale_fill_manual(  values = arm_colors) +
    scale_y_reverse() +
    labs(title = "A.  Disease Score Over Time",
         subtitle = "Spaghetti + mean \u00b1 SEM (Mixed-Effects Model)",
         x = "Week", y = "Disease Score (0\u201310, lower = better)",
         colour = NULL, fill = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── B: IgG PD response with Emax curve fits ───────────────
  active_long <- longitudinal |>
    dplyr::filter(treatment == 1) |>
    dplyr::mutate(resp_label = ifelse(true_responder == 1,
                                      "Responders", "Non-Responders"))

  igg_mean <- active_long |>
    dplyr::group_by(resp_label, week) |>
    dplyr::summarise(
      mean_pct = mean(pct_igg_change),
      sem_pct  = sd(pct_igg_change) / sqrt(dplyr::n()),
      .groups  = "drop"
    )

  emax_fn <- function(w, emax, ec50) -emax * w / (ec50 + w)
  fit_weeks <- seq(0, 24, by = 0.5)
  curve_df  <- purrr::map_dfr(names(emax_res$params), function(nm) {
    p <- emax_res$params[[nm]]
    tibble::tibble(
      resp_label = nm,
      week       = fit_weeks,
      fit_pct    = emax_fn(fit_weeks, p["emax"], p["ec50"])
    )
  })

  resp_colors <- c("Responders" = "#d62728", "Non-Responders" = "#ff7f0e")

  pB <- ggplot(igg_mean, aes(week, mean_pct, colour = resp_label)) +
    geom_ribbon(aes(ymin = mean_pct - sem_pct,
                    ymax = mean_pct + sem_pct, fill = resp_label),
                alpha = 0.15, colour = NA) +
    geom_line(data = curve_df, aes(y = fit_pct),
              linewidth = 1.5, linetype = "dashed") +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = resp_colors) +
    scale_fill_manual(  values = resp_colors) +
    labs(title = "B.  IgG PD Response — Active Arm",
         subtitle = "Emax curve fits by responder status",
         x = "Week", y = "% Change from Baseline IgG",
         colour = NULL, fill = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── C: Kaplan-Meier curves ────────────────────────────────
  km_df <- with(
    survfit(Surv(time_to_response_weeks, event_observed) ~ treatment,
            data = survival_df),
    tibble::tibble(
      time   = time,
      surv   = surv,
      lower  = lower,
      upper  = upper,
      arm    = rep(names(strata), strata)
    )
  ) |>
    dplyr::mutate(arm = ifelse(grepl("treatment=1", arm), "Active", "Placebo"))

  # Add t=0 start
  km_start <- tibble::tibble(
    time = 0, surv = 1, lower = 1, upper = 1,
    arm = c("Active","Placebo")
  )
  km_plot_df <- dplyr::bind_rows(km_start, km_df)

  p_label <- sprintf("Log-rank p = %.5f", km_res$lr_p)

  pC <- ggplot(km_plot_df, aes(time, surv, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_step(linewidth = 1.6) +
    scale_colour_manual(values = arm_colors) +
    scale_fill_manual(  values = arm_colors) +
    annotate("text", x = 15, y = 0.90, label = p_label, size = 3.5) +
    labs(title = "C.  Kaplan-Meier: Time-to-Response",
         x = "Weeks", y = "P(no response yet)",
         colour = NULL, fill = NULL) +
    coord_cartesian(ylim = c(0, 1.05)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── D: Cox forest plot ────────────────────────────────────
  forest_df <- cox_res$hr_df |>
    dplyr::mutate(
      term = dplyr::recode(term,
        treatment        = "Treatment (Active vs Placebo)",
        baseline_igg_z   = "Baseline IgG (standardised)",
        latent_biology_z = "Latent biology score (standardised)"
      ),
      direction = ifelse(HR >= 1, "Faster response", "Slower response")
    )

  pD <- ggplot(forest_df, aes(y = reorder(term, HR))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper,
                       colour = direction), height = 0.25, linewidth = 1.2) +
    geom_point(aes(x = HR, colour = direction), size = 4) +
    scale_x_log10() +
    scale_colour_manual(values = c("Faster response" = "#d62728",
                                   "Slower response" = "#1f77b4")) +
    labs(title = "D.  Cox PH Forest Plot",
         subtitle = "HR > 1 = faster time-to-response",
         x = "Hazard Ratio (95% CI)  [log scale]",
         y = NULL, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  # ── Combine and save ──────────────────────────────────────
  combined <- (pA + pB) / (pC + pD) +
    patchwork::plot_annotation(
      title = "Module 2 — Longitudinal Modelling & Survival Analysis",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )

  out_path <- file.path(output_dir, "02_longitudinal_survival.png")
  ggsave(out_path, combined, width = 14, height = 11, dpi = 200)
  cat(sprintf("\n  \u2713 Saved: %s\n", out_path))
}
