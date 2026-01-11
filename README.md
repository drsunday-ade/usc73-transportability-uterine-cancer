# USC73 Transportability in Uterine Serous Carcinoma (TCGA ↔ UPSC)

**Repository:** `usc73-transportability-uterine-cancer`  
**Focus:** Transportability (bi-directional internal–external validation) of **fixed-horizon survival prediction** in uterine serous carcinoma under **cohort shift**, comparing a clinicopathologic model **with vs without** the **USC73 molecular signature**.

> Research code for statistical evaluation of prediction models. Not for clinical decision-making.

---

## 1) What this project answers

**Primary study question:**  
When a fixed-horizon survival prediction model is developed in one uterine serous carcinoma cohort and transported to another, **how well does it generalize**, and **does adding USC73 provide incremental, transportable value** beyond standard clinicopathologic predictors?

**Design feature (core contribution):**  
A **bi-directional internal–external validation** (TCGA→UPSC and UPSC→TCGA) to ensure conclusions reflect **true out-of-cohort performance** under shift (not within-cohort optimism).

---

## 2) Key results (from the pipeline outputs)

- Total **N (TCGA + UPSC) = 125**; complete-case age **N = 122**.  
- Candidate horizons: **τ = 1–5 years**; representative horizon **τ\* ≈ 2 years** (both transport directions).  
- Best model by **mean recalibrated Brier**:
  - **TCGA→UPSC:** *Cox Base* — mean recalibrated Brier **0.1394**, mean C-index **0.7714**.
  - **UPSC→TCGA:** *Ridge Base* — mean recalibrated Brier **0.1860**, mean C-index **0.7495**.
- Incremental value (Extended − Base) at τ\* shows **directional asymmetry**:
  - **TCGA→UPSC (Cox):** ΔC-index **−0.093**; ΔBrier\_recal **+0.0240** (worse).
  - **UPSC→TCGA (Cox):** ΔC-index **+0.044**; ΔBrier\_recal **−0.0131** (improves).
- **Time-dependent AUC** can be **non-computable** at some τ due to estimability/weight-stability safeguards; such settings are **flagged rather than forced** into summaries.

See: `results/results_summary.txt` and the **Key** tables under `tables/`.

---

## 3) Repository structure (run artifacts)

This repository is organized around a single, fully scripted pipeline run directory structure:

```
.
├── data/                # input data (not distributed publicly; see Data section)
├── scripts/             # analysis pipeline(s)
├── figures/             # publication figures (PNG)
├── tables/              # publication tables (CSV + HTML; plus _extras)
├── results/             # compact summaries (e.g., results_summary.txt)
├── overview-of-data/    # variable dictionaries / manifests produced by the pipeline
└── logs/                # pipeline_log.txt, sessionInfo.txt, manifests
```

---

## 4) Input data expectations

The pipeline reads an Excel workbook (default: `Supplementary_Expression_phenotype_data.xlsx`) and either:
- **auto-detects** the appropriate sheet and key columns, or
- uses `--sheet` to specify the sheet explicitly.

### Required information (minimum)
The code expects patient-level data with:
- `sample_id` (used to infer cohort membership via `TCGA` vs `UPSC`)
- survival **time** and **event** indicator
- clinicopathologic covariates (see below)
- `usc73_score` for extended-model analyses

The script uses **pattern-based column detection** (robust to minor naming differences). If detection fails, rename columns to include common substrings such as `time`, `os`, `status`, `event`, `age`, `stage`, `race`, `usc73`.

---

## 5) Covariates used (extracted from the analysis script)

**Base model predictors (prespecified):**
- `age_at_dx` (continuous)
- `race3` (3-level factor: **AA**, **C**, **Other/Unknown**)  
- `stage` (binary: **early** vs **advanced**)

**Extended model predictors:**
- all base predictors **plus**
- `usc73_score` (continuous)

These definitions come directly from `scripts/USC73_transportability_pipeline_v8.R` (`base_vars` and `ext_vars`).

---

## 6) Methods implemented (high-level)

Within each transport direction (TCGA→UPSC and UPSC→TCGA), the pipeline evaluates:

- **Model families:** Cox PH and ridge-penalized Cox (glmnet)
- **Predictor sets:** Base vs Extended (USC73 added)
- **Transport vs recalibration:** raw transported predictions vs **cross-fitted** target recalibration
- **Censoring-aware metrics:** IPCW Brier score, C-index, calibration (slope/intercept), and **IPCW decision curve analysis**
- **Estimability gating:** prevents reporting unstable/non-identifiable fixed-horizon metrics when follow-up support is inadequate

---

## 7) Quickstart

### A) Requirements
- R (recommended: R ≥ 4.2)
- Core packages are installed automatically by default (toggle with `--no_install TRUE`)

Core packages used:
- `survival`, `glmnet`
- `ggplot2`, `dplyr`, `tidyr`, `tibble`
- `readxl`, `writexl`
- `htmlTable`
Optional (for enhanced plots/metrics): `viridisLite`, `ggrepel`, `timeROC`

### B) Run the pipeline
From the repository root:

```bash
Rscript scripts/USC73_transportability_pipeline_v8.R \
  --input data/Supplementary_Expression_phenotype_data.xlsx \
  --outdir usc73_v8_run
```

Useful options:

```bash
# choose the Excel sheet explicitly
--sheet "Sheet1"

# reproducibility
--seed 2025

# cross-fitted recalibration folds
--nfolds 5

# candidate horizons (years)
--tau_candidates "1,2,3,4,5"

# disable automatic package installation
--no_install TRUE

# write HTML versions of the Key tables
--make_html TRUE
```

### C) Outputs you should expect
The pipeline writes:
- **Figures:** `figures/Fig01_...png` to `figures/Fig18_...png`
- **Tables:** `tables/Table_K0_...` to `tables/Table_K12_...` (CSV + optional HTML)
- **Logs:** `logs/pipeline_log.txt`, `logs/sessionInfo.txt`
- **Summary:** `results/results_summary.txt`

---

## 8) Outputs: what to look at first

If you want a reviewer-friendly pass through the evidence:

1. **Cohort shift & outcome context**  
   - `Fig01_CohortShift_Bars.png`, `Fig03_KM_ByCohort.png`, `Fig15_StageDistribution_ByCohort.png`, `Fig16_AgeDistribution_ByCohort.png`

2. **Transport performance across horizons (τ = 1–5 years)**  
   - `Fig07_CIndex_vs_Tau.png`, `Fig08_Brier_vs_Tau_Raw_vs_Recal.png`, `Fig09_DeltaCIndex_ExtMinusBase.png`, `Fig14_CalibrationParams_vs_Tau.png`

3. **Representative horizon (τ* ≈ 2 years): calibration + utility**  
   - `Fig10_CalibrationCurve_RepTau.png`, `Fig11_DecisionCurves_RepTau.png`, `Fig18_DCA_Gain_vs_All.png`, `Fig12_RiskDistribution_RepTau.png`

4. **Key tables for the manuscript**  
   - `Table_K6_Performance_Key.*` (core performance summaries)
   - `Table_K7_Deltas_ExtMinusBase_Key.*` (incremental value, Extended − Base)
   - `Table_K5_Estimability_Key.*` (what was estimable and why)

---

## 9) Reproducibility and quality guardrails

This pipeline is designed to be **defensive** under cohort shift and censoring:

- **Strict cohort separation** for development vs evaluation
- **No stepwise selection** (prespecified covariates)
- **IPCW** for all fixed-horizon evaluation under right censoring
- **Estimability gating**: flags non-computable settings rather than producing misleading metrics
- **Cross-fitted recalibration** in the target cohort to reduce optimism when learning calibration mappings
- **Session capture** via `logs/sessionInfo.txt`

---

## 10) Data availability

The original cohort data sources may have access restrictions. This repository does **not** distribute individual-level patient data by default. To reproduce:
1. Obtain the Excel workbook used in your analysis environment.
2. Place it under `data/` (or point `--input` to its location).
3. Run the pipeline as shown above.

If you plan to share a de-identified or synthetic example dataset, add it to `data/` and document it here.

---

## 11) How to cite

If you use this code in academic work, please cite:
- The associated manuscript (add citation here once submitted/accepted)
- This repository (GitHub) and/or a Zenodo DOI (recommended)

**Recommended additions**
- `CITATION.cff` (GitHub will render this automatically)
- Zenodo archiving for a DOI on release

---

## 12) License and disclaimer

- **License:** add your license of choice (MIT/Apache-2.0 are common for research code).
- **Disclaimer:** This repository is for **research and educational** use only and is not validated for clinical deployment.

---

## 13) Contact

Maintainer: **Sunday Adetunji**  
For issues, please open a GitHub Issue with:
- your `logs/pipeline_log.txt`
- your `logs/sessionInfo.txt`
- the exact command used
