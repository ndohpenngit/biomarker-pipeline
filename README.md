# Biomarker Data Science Platform

End-to-end clinical trial biomarker analysis pipelines across four therapeutic areas, built in R and published as a Quarto website.

**Live site →** `https://ndohpenngit.github.io/biomarker-pipeline/`

---

## Therapeutic Areas

| Report | Disease | Mechanism | Primary Biomarker | N |
|---|---|---|---|---|
| [Immunology](biomarker_reports/autoimmune_fcrn.qmd) | Generalised Myasthenia Gravis | FcRn Inhibitor (VYVGART-like) | IgG total NPX | 120 |
| [Oncology](biomarker_reports/oncology_checkpoint.qmd) | NSCLC | PD-1 Checkpoint Inhibitor (nivolumab-like) | TMB | 150 |
| [Cardiovascular](biomarker_reports/cardiovascular_pcsk9.qmd) | Heterozygous FH | PCSK9 Inhibitor (evolocumab-like) | LDL-C | 200 |
| [Neurology](biomarker_reports/neurology_alzheimers.qmd) | Early Alzheimer's Disease | Anti-Aβ mAb (lecanemab-like) | p-tau181 | 160 |

All data are **fully simulated** (fixed seeds per report). No patient data are used.

---

## Pipeline

Each report runs an identical three-stage analytical pipeline parameterised per therapeutic area:

```
_analysis_core.qmd          ← shared source (included, not rendered directly)
├── Step 1 · Multi-Omics
│     Transcriptomics QC · batch correction · Welch DE (BH-FDR)
│     Olink NPX proteomics · cross-modal correlation
├── Step 2 · Longitudinal & Survival
│     lme4 linear mixed-effects · Emax PD model (nls)
│     Kaplan–Meier · Cox proportional-hazards
└── Step 3 · Machine Learning
      PCA + UMAP · k-means clustering
      Elastic-net (glmnet) · Random forest · OOB AUROC
```

TA-specific framing (disease background, biomarker tables, clinical context) lives in each wrapper `.qmd`; all R code lives in `_analysis_core.qmd`.

---

## Project Structure

```
biomarker_reports/
├── _quarto.yml                  # Quarto website project config
├── _analysis_core.qmd           # Shared parameterised analysis engine
├── index.qmd                    # Landing page
├── autoimmune_fcrn.qmd          # Immunology report
├── oncology_checkpoint.qmd      # Oncology report
├── cardiovascular_pcsk9.qmd     # Cardiovascular report
├── neurology_alzheimers.qmd     # Neurology report
├── render_all.R                 # Render helper
├── .gitignore
├── docs/                        # Rendered site (GitHub Pages source)
└── R/
    ├── 01_simulate_data.R
    ├── 02_multiomics_analysis.R
    ├── 03_longitudinal_survival.R
    └── 04_ml_pipeline.R
```

---

## Quickstart

### 1 — Install dependencies

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "patchwork",
  "lme4", "lmerTest", "survival", "broom", "broom.mixed",
  "glmnet", "randomForest", "cluster", "pROC", "umap",
  "knitr", "kableExtra", "sessioninfo",
  "survminer", "forcats", "ggrepel", "pheatmap"
))

# Bioconductor
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("limma", "SummarizedExperiment"))
```

### 2 — Preview locally

```bash
cd biomarker_reports
quarto preview          # opens site at localhost in browser
```

### 3 — Render and publish

```bash
quarto render           # builds to docs/
git add docs/
git commit -m "Render site"
git push                # GitHub Pages updates automatically
```

---

## Tech Stack

| Layer | Tool |
|---|---|
| Reporting | [Quarto](https://quarto.org) |
| Language | R 4.x |
| Visualisation | ggplot2, patchwork, pheatmap |
| Statistics | lme4, survival, limma |
| Machine learning | glmnet, randomForest, pROC, umap |
| Publishing | GitHub Pages |
