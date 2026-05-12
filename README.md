# Biomarker Data Science Pipeline (R)
### End-to-End Autoimmune Trial Biomarker Analysis · FcRn Inhibitor (VYVGART-like)

---

## Overview

A complete R-based biomarker analysis pipeline demonstrating the statistical
and computational skills required for immunology drug development. The pipeline
simulates a Phase 2 randomised trial of an FcRn inhibitor in a severe autoimmune
disease (e.g. generalised Myasthenia Gravis, gMG), and covers multi-omics QC,
longitudinal modelling, survival analysis, and predictive machine learning —
all in idiomatic, tidyverse-fluent R.

```
biomarker_pipeline_r/
├── run_pipeline.R              # Single entry point
├── R/
│   ├── 01_simulate_data.R      # Trial data generator
│   ├── 02_multiomics_analysis.R # Module 1: transcriptomics + proteomics
│   ├── 03_longitudinal_survival.R # Module 2: lme4, Emax, KM, Cox PH
│   └── 04_ml_pipeline.R        # Module 3: PCA, K-means, glmnet, randomForest
├── outputs/
│   ├── 01_multiomics_analysis.png
│   ├── 02_longitudinal_survival.png
│   └── 03_ml_pipeline.png
└── README.md
```

---

## Simulated Biology

| Signal | Implementation |
|--------|---------------|
| FcRn-mediated IgG catabolism | Emax-style reduction in active arm |
| Responder phenotype | ~70% active / ~25% placebo response rate |
| Transcriptomic signature | 50 baseline DE genes up-regulated in responders |
| Batch effect | 2 sequencing runs (median-centering correction) |
| Proteomics (Olink NPX) | 50 proteins × 3 timepoints; IgG1/4/total reduced by FcRn blockade |
| Time-to-response | Weibull-distributed; active arm responds earlier |
| Censoring | ~10% dropout |

---

## Module 1 — Multi-Omics (`R/02_multiomics_analysis.R`)

### Transcriptomics QC
1. **log₂(count + 1)** normalisation
2. **Median-centering batch correction** (mirrors `limma::removeBatchEffect`
   for a 2-batch design — loadable when Bioconductor is available)
3. **High-variance gene filter** — top 50% by per-gene variance
4. **PCA outlier detection** — samples > 97.5th percentile Euclidean distance

### Differential Expression
- Per-gene **Welch t-test** + **Benjamini-Hochberg FDR** (base R, no external dependency)
- Designed to swap in `limma::eBayes` or `DESeq2` without changing the interface
- Results: 30 significant baseline DE genes (FDR < 5%)

### Proteomics (Olink NPX)
- Longitudinal NPX Δ (Wk12 – Wk0) per protein per arm
- IgG_total: Active –1.79 NPX vs Placebo –0.09 NPX

### Multi-Omics Integration
- Transcriptomic PCA → Pearson r vs key baseline proteins
- Production path: `MOFA+`, `mixOmics::DIABLO`, or weighted SNF

---

## Module 2 — Longitudinal & Survival (`R/03_longitudinal_survival.R`)

### Linear Mixed-Effects Model (`lme4`)
```r
lmer(disease_score ~ week * treatment + (1 + week | patient_id))
```
- **week:treatment** interaction = primary drug-effect estimand
- Result: –1.45 disease-score units drug effect at Week 24
- `lmerTest` adds Satterthwaite df and p-values when available

### Emax Pharmacodynamic Model (`nls`)
```
ΔIgG% = –Emax × week / (EC50 + week)
```
- Responders: **Emax = 83%, EC50 = 4.1 weeks**
- Non-Responders: **Emax = 37%, EC50 = 4.0 weeks**

### Kaplan-Meier + Log-rank (`survival`)
- Log-rank χ² = 30.1, **p < 0.0001**
- Active median time-to-response: 14.3 weeks

### Cox Proportional Hazards
| Covariate | HR | Interpretation |
|-----------|-----|----------------|
| treatment | 10.5× | Active arm responds ~10× faster |
| baseline_igg_z | 1.0× | Minimal effect |
| latent_biology_z | 2.3× | Biology strongly predicts response |

---

## Module 3 — ML Pipeline (`R/04_ml_pipeline.R`)

### Dimensionality Reduction
- **PCA** (centred + scaled) → n PCs for 90% variance
- **UMAP** (`umap` package, n_neighbors=15, min_dist=0.1) on top 20 PCs
  - Falls back to t-SNE (`Rtsne`) or PC1/PC2 if `umap` unavailable

### Patient Clustering
- **K-means** in top-10 PC space; K via elbow analysis
- Silhouette score for internal validation
- Cluster ↔ responder rate as post-hoc validation

### Predictive Biomarker Modelling (Two-Stage)
**Stage 1 — `glmnet` ElasticNet** (multi-alpha CV grid):
- Combined matrix: 250 high-variance genes + 50 proteins = 300 features
- Sparse selection → biologically interpretable panel

**Stage 2 — `randomForest`** (5-fold stratified CV):
- Class-weighted for responder imbalance
- **CV AUROC = 0.993 ± 0.016**
- OOB AUROC = 0.979
- Feature importance: Mean Decrease in Gini

---

## Key Findings (Simulated Data)

```
Multi-Omics:
  • 30 baseline DE genes (Welch t-test, FDR < 5%)
  • IgG_total Δ: Active –1.79 NPX vs Placebo –0.09 NPX (Wk0→12)

Longitudinal & Survival:
  • lme4 week×treatment: β = –0.06/wk  → –1.45 units at Week 24
  • Emax: Responders 83% vs Non-Responders 37% IgG reduction
  • Log-rank p < 0.0001 — active arm responds earlier
  • Cox treatment HR = 10.5×

ML:
  • glmnet selected features from 300-feature matrix
  • randomForest OOB AUROC = 0.979
```

---

## Running the Pipeline

```r
# Install dependencies (one time)
install.packages(c("ggplot2","dplyr","tidyr","purrr","tibble",
                   "lme4","lmerTest","survival","broom",
                   "glmnet","randomForest","cluster",
                   "patchwork","ggrepel","pheatmap","RColorBrewer",
                   "pROC","umap"))

# Run full pipeline (~60–90 seconds)
Rscript run_pipeline.R
```

---

## R Package Stack vs argenx JD

| JD Requirement | R Implementation |
|----------------|-----------------|
| Longitudinal models | `lme4::lmer` + `lmerTest` |
| Survival analysis | `survival::coxph` + KM |
| Multiplicity control | BH FDR (base `p.adjust`) |
| High-dimensional DE | Welch t-test; swap to `limma`/`DESeq2` |
| Dimensionality reduction | `prcomp` + `umap`/`Rtsne` |
| Clustering | `stats::kmeans` + `cluster::silhouette` |
| Regularised models | `glmnet` (ElasticNet) |
| ML for biomarker discovery | `randomForest` + `pROC` |
| Reproducible code | `set.seed`, modular R functions |
| Visualisation | `ggplot2` + `patchwork` + `pheatmap` |
| Dashboarding (prod) | R/Shiny, Posit Connect |
| Cloud platforms (prod) | Databricks, Snowflake, Posit Team |
