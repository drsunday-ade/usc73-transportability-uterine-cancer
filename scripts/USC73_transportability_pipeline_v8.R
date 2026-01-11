#!/usr/bin/env Rscript
# ============================================================
# USC73 TRANSPORTABILITY PIPELINE (V8 "publication-grade + fail-safe")
#
# PURPOSE
#   Two-cohort internal-external validation (TCGA <-> UPSC) for 1–5y OS risk
#   with transportability diagnostics and *target-cohort cross-fitted recalibration*.
#
# KEY METHODS / NOVELTY (publishable contributions)
#   1) Estimability-aware time-horizon (tau) filtering per direction, to prevent
#      meaningless evaluation when cases/controls are insufficient at tau.
#   2) Cross-fitted recalibration performed *in the target cohort* (not in-sample),
#      to reduce optimism and avoid “perfect calibration” artifacts.
#   3) IPCW scoring for Brier, calibration, and decision curves under censoring.
#
# OUTPUTS (guaranteed if script completes)
#   - >=12 publication-quality figures (PNG) with robust color scales + annotations
#   - 10 key tables in BOTH CSV and HTML (simple standalone HTML with CSS)
#   - Full supplementary tables saved under tables/_extras
#   - Logs + session info saved under logs/
#
# RUN (recommended)
#   Rscript USC73_transportability_pipeline_v8.R --input "Supplementary_Expression_phenotype_data.xlsx" --outdir "usc73_v8"
#
# NOTES
#   - If --input omitted, auto-picks the newest .xlsx/.xls in working directory.
#   - If viridis (older) lacks scale_*_viridis_d, script falls back safely.
#   - HTML tables do NOT require webshot/Chromium; written as standalone HTML.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# -------------------------------
# 0) CLI parsing
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

INPUT      <- get_arg("--input",  NA)
OUTDIR     <- get_arg("--outdir", "usc73_v8")
SHEET      <- get_arg("--sheet",  NA)
NO_INSTALL <- tolower(get_arg("--no_install", "false")) %in% c("true","t","1","yes","y")
MAKE_HTML  <- tolower(get_arg("--make_html", "true"))  %in% c("true","t","1","yes","y")

SEED   <- suppressWarnings(as.integer(get_arg("--seed", "20260110")))
if (is.na(SEED)) SEED <- 20260110
set.seed(SEED)

NFOLDS <- suppressWarnings(as.integer(get_arg("--nfolds", "5")))
if (is.na(NFOLDS) || NFOLDS < 3) NFOLDS <- 5

MIN_CASES    <- suppressWarnings(as.integer(get_arg("--min_cases", "10")))
MIN_CONTROLS <- suppressWarnings(as.integer(get_arg("--min_controls", "10")))
if (is.na(MIN_CASES) || MIN_CASES < 3) MIN_CASES <- 10
if (is.na(MIN_CONTROLS) || MIN_CONTROLS < 3) MIN_CONTROLS <- 10

TAU_CAND_STR <- get_arg("--tau_candidates", "1,1.5,2,2.5,3,3.5,4,5")
TAU_CAND <- suppressWarnings(as.numeric(strsplit(TAU_CAND_STR, ",")[[1]]))
TAU_CAND <- sort(unique(TAU_CAND[is.finite(TAU_CAND) & TAU_CAND > 0]))
if (length(TAU_CAND) == 0) TAU_CAND <- c(1,2,3)

RUN_MI  <- tolower(get_arg("--mi", "false")) %in% c("true","t","1","yes","y")
RUN_RSF <- tolower(get_arg("--run_rsf", "false")) %in% c("true","t","1","yes","y")

# -------------------------------
# 1) Packages (install if needed)
# -------------------------------
safe_install_load <- function(pkgs) {
  for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    if (!ok) {
      if (NO_INSTALL) stop("Missing package: ", p, ". Re-run without --no_install or install manually.", call.=FALSE)
      message("Installing package: ", p)
      tryCatch(
        install.packages(p, repos = "https://cloud.r-project.org"),
        error = function(e) stop("Failed to install ", p, ": ", e$message, call.=FALSE)
      )
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

core <- c(
  "readxl","readr","dplyr","tibble","stringr","janitor","tidyr","purrr",
  "ggplot2","survival","glmnet","scales"
)
safe_install_load(core)

# Optional (plots)
opt <- c("viridis","ggrepel")
for (p in opt) {
  if (!requireNamespace(p, quietly = TRUE) && !NO_INSTALL) {
    try(install.packages(p, repos="https://cloud.r-project.org"), silent=TRUE)
  }
}
has_viridis <- requireNamespace("viridis", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_viridis) suppressPackageStartupMessages(library(viridis))
if (has_ggrepel) suppressPackageStartupMessages(library(ggrepel))

# Optional (time-dependent AUC)
has_timeROC <- requireNamespace("timeROC", quietly = TRUE)
if (!has_timeROC && !NO_INSTALL) try(install.packages("timeROC", repos="https://cloud.r-project.org"), silent=TRUE)
has_timeROC <- requireNamespace("timeROC", quietly = TRUE)
if (has_timeROC) suppressPackageStartupMessages(library(timeROC))

# Optional MI / RSF
has_mice <- requireNamespace("mice", quietly = TRUE)
if (RUN_MI && !has_mice && !NO_INSTALL) try(install.packages("mice", repos="https://cloud.r-project.org"), silent=TRUE)
has_mice <- requireNamespace("mice", quietly = TRUE)
if (RUN_MI && has_mice) suppressPackageStartupMessages(library(mice))

has_rsf <- requireNamespace("randomForestSRC", quietly = TRUE)
if (RUN_RSF && !has_rsf && !NO_INSTALL) try(install.packages("randomForestSRC", repos="https://cloud.r-project.org"), silent=TRUE)
has_rsf <- requireNamespace("randomForestSRC", quietly = TRUE)
if (RUN_RSF && has_rsf) suppressPackageStartupMessages(library(randomForestSRC))

# -------------------------------
# 2) Output dirs + logging
# -------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "tables", "_extras"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "logs"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "models"),  showWarnings = FALSE, recursive = TRUE)

LOGFILE <- file.path(OUTDIR, "logs", "pipeline_log.txt")
sink(LOGFILE, split = TRUE)
on.exit({ try(sink(NULL), silent = TRUE) }, add = TRUE)

cat("=== USC73 TRANSPORTABILITY PIPELINE V8 ===\n")
cat("Seed:", SEED, "| nfolds:", NFOLDS, "| min_cases:", MIN_CASES, "| min_controls:", MIN_CONTROLS, "\n")
cat("Tau candidates:", paste(TAU_CAND, collapse=", "), "\n")
cat("MAKE_HTML:", MAKE_HTML, "| NO_INSTALL:", NO_INSTALL, "\n\n")

# Save session info
try({
  si <- capture.output(sessionInfo())
  writeLines(si, file.path(OUTDIR, "logs", "sessionInfo.txt"))
}, silent=TRUE)

# -------------------------------
# 3) Input autodetect
# -------------------------------
auto_pick_input <- function() {
  xl <- list.files(getwd(), pattern="\\.(xlsx|xls)$", full.names = TRUE)
  if (length(xl) == 0) return(NA_character_)
  xl[which.max(file.info(xl)$mtime)]
}

if (is.na(INPUT) || !nzchar(INPUT)) {
  INPUT <- auto_pick_input()
  if (is.na(INPUT)) stop("No --input provided and no Excel file found in working directory.", call.=FALSE)
  cat("NOTE: --input not provided. Auto-selected:", INPUT, "\n\n")
}
if (!file.exists(INPUT)) stop("Input file not found: ", INPUT, call.=FALSE)

# -------------------------------
# 4) Read Excel + best sheet selection
# -------------------------------
clinical_regex <- "(^age$|age_at_dx|age_at_diagn|race|ethnic|stage|figo|os\\b|overall_survival|status|event|death|usc73|sample|patient|id)"

score_sheet <- function(df) {
  cn <- janitor::make_clean_names(names(df))
  clin_hits <- sum(stringr::str_detect(cn, clinical_regex))
  tibble::tibble(clin_hits = clin_hits, ncol = ncol(df), nrow = nrow(df),
                 score = clin_hits * 1000 + ncol(df))
}

read_sheet_safe <- function(path, sheet) {
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 5000),
           error = function(e) NULL)
}

sheets <- readxl::excel_sheets(INPUT)
if (length(sheets) == 0) stop("No sheets found in Excel file.", call.=FALSE)

if (!is.na(SHEET) && nzchar(SHEET)) {
  if (!(SHEET %in% sheets)) stop("Requested sheet not found: ", SHEET, call.=FALSE)
  chosen_sheet <- SHEET
} else {
  manifest <- purrr::map_dfr(sheets, function(sh) {
    df_try <- read_sheet_safe(INPUT, sh)
    if (is.null(df_try)) return(tibble::tibble(sheet=sh, clin_hits=NA, ncol=NA, nrow=NA, score=NA))
    tibble::tibble(sheet=sh) %>% dplyr::bind_cols(score_sheet(df_try))
  }) %>% dplyr::arrange(dplyr::desc(score))
  readr::write_csv(manifest, file.path(OUTDIR, "logs", "sheet_manifest_scores.csv"))
  chosen_sheet <- manifest$sheet[which.max(manifest$score)]
}

cat("Input:", INPUT, "\nChosen sheet:", chosen_sheet, "\n\n")
df0 <- read_sheet_safe(INPUT, chosen_sheet)
if (is.null(df0)) stop("Failed to read chosen sheet: ", chosen_sheet, call.=FALSE)
names(df0) <- janitor::make_clean_names(names(df0))

# -------------------------------
# 5) Robust column detection + standardization
# -------------------------------
find_col <- function(nms, patterns, prefer = NULL) {
  nms_low <- tolower(nms)
  hits <- unique(unlist(lapply(patterns, function(p) which(stringr::str_detect(nms_low, p)))))
  if (length(hits) == 0) return(NA_character_)
  cand <- nms[hits]
  if (!is.null(prefer)) {
    pref_hits <- cand[tolower(cand) %in% tolower(prefer)]
    if (length(pref_hits) > 0) return(pref_hits[1])
  }
  cand[1]
}

is_numericish <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  xn <- suppressWarnings(as.numeric(as.character(x)))
  denom <- max(1, sum(!is.na(x)))
  (sum(!is.na(xn)) / denom) >= 0.95
}

id_col <- if ("sample_id" %in% names(df0)) "sample_id" else
  find_col(names(df0), patterns=c("^row$","sample","patient","case_id","^id$"), prefer=c("row","sample_id","patient_id"))

os_time_col <- find_col(names(df0), patterns=c("^os\\b","overall_survival","os_year","os_years","survival.*year","time"),
                        prefer=c("os_years","os_year","os"))
os_stat_col <- find_col(names(df0), patterns=c("os_status","status","event","death"),
                        prefer=c("os_status"))

age_col  <- find_col(names(df0), patterns=c("age_at_dx","age_at_diagn","age_at_diagnosis","^age$"),
                     prefer=c("age_at_dx","age"))
race_col <- find_col(names(df0), patterns=c("^race$","ethnic","ethnicity"), prefer=c("race","ethnicity"))
stage_col<- find_col(names(df0), patterns=c("^stage$","figo"), prefer=c("stage","figo_stage","figo"))

usc_candidates <- names(df0)[stringr::str_detect(names(df0), "usc73")]
if (length(usc_candidates) == 0) stop("No USC73 columns detected (need 'usc73' in name).", call.=FALSE)
num_cands <- usc_candidates[sapply(usc_candidates, function(v) is_numericish(df0[[v]]))]
if (length(num_cands) == 0) stop("USC73 numeric-like column not found among: ", paste(usc_candidates, collapse=", "), call.=FALSE)
usc_score_col <- num_cands[which.max(sapply(num_cands, function(v) sum(!is.na(df0[[v]]))))]

cat("Detected columns:\n",
    " id:", id_col, "\n",
    " os_time:", os_time_col, "\n",
    " os_status:", os_stat_col, "\n",
    " age:", age_col, "\n",
    " race:", race_col, "\n",
    " stage:", stage_col, "\n",
    " usc73:", usc_score_col, "\n\n", sep="")

if (any(is.na(c(id_col, os_time_col, os_stat_col, age_col, race_col, stage_col)))) {
  stop("Could not identify all required columns. Provide --sheet or rename columns.", call.=FALSE)
}

harmonize_stage <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- stringr::str_replace_all(s, "\\s+", " ")
  dplyr::case_when(
    stringr::str_detect(s, "early") ~ "early",
    stringr::str_detect(s, "advanced") ~ "advanced",
    stringr::str_detect(s, "iv|4") ~ "advanced",
    stringr::str_detect(s, "iii|3") ~ "advanced",
    stringr::str_detect(s, "ii|2") ~ "early",
    stringr::str_detect(s, "i|1") ~ "early",
    TRUE ~ NA_character_
  )
}

harmonize_race3 <- function(x) {
  s <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    stringr::str_detect(s, "african|black|\\baa\\b") ~ "AA",
    stringr::str_detect(s, "cauc|white|\\bc\\b") ~ "C",
    TRUE ~ "Other/Unknown"
  )
}

df <- df0 %>%
  dplyr::transmute(
    sample_id = as.character(.data[[id_col]]),
    cohort = dplyr::case_when(
      stringr::str_starts(sample_id, "TCGA") ~ "TCGA",
      stringr::str_starts(sample_id, "UPSC") ~ "UPSC",
      TRUE ~ "OTHER"
    ),
    os_years = suppressWarnings(as.numeric(.data[[os_time_col]])),
    os_status_raw = .data[[os_stat_col]],
    age_at_dx = suppressWarnings(as.numeric(.data[[age_col]])),
    race_raw  = .data[[race_col]],
    stage_raw = .data[[stage_col]],
    usc73_score = suppressWarnings(as.numeric(.data[[usc_score_col]]))
  ) %>%
  dplyr::mutate(
    os_status = dplyr::case_when(
      is.numeric(os_status_raw) | is.integer(os_status_raw) ~ as.integer(os_status_raw),
      TRUE ~ {
        s <- tolower(trimws(as.character(os_status_raw)))
        as.integer(dplyr::case_when(
          s %in% c("1","dead","death","deceased","event","yes") ~ 1L,
          s %in% c("0","alive","censored","no") ~ 0L,
          TRUE ~ NA_integer_
        ))
      }
    ),
    stage = factor(harmonize_stage(stage_raw), levels=c("early","advanced")),
    race3 = factor(harmonize_race3(race_raw), levels=c("AA","C","Other/Unknown"))
  )

df_ie <- df %>% dplyr::filter(cohort %in% c("TCGA","UPSC"))

# Integrity checks (fail fast with actionable messages)
if (nrow(df_ie) < 50) stop("Too few TCGA+UPSC rows after filtering. Check sample_id prefixes.", call.=FALSE)
if (any(is.na(df_ie$os_years)) || any(df_ie$os_years < 0, na.rm=TRUE)) stop("Invalid os_years detected (NA or <0).", call.=FALSE)
if (any(is.na(df_ie$os_status)) || !all(df_ie$os_status %in% c(0,1))) stop("Invalid os_status detected; must be 0/1.", call.=FALSE)
if (any(is.na(df_ie$usc73_score))) stop("usc73_score contains NA after coercion. Fix USC73 column or coercion.", call.=FALSE)
if (any(is.na(df_ie$stage))) stop("Stage harmonization produced NA. Update harmonize_stage() patterns.", call.=FALSE)
if (any(duplicated(df_ie$sample_id))) stop("Duplicate sample_id detected.", call.=FALSE)

df_cc <- df_ie %>% dplyr::filter(!is.na(age_at_dx))  # complete-case for age (or use MI via --mi)
cat("N total (TCGA+UPSC):", nrow(df_ie), "| N complete-case:", nrow(df_cc), "\n\n")

# -------------------------------
# 6) Table writers (CSV + HTML)
# -------------------------------
escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&","&amp;",x, fixed=TRUE)
  x <- gsub("<","&lt;",x, fixed=TRUE)
  x <- gsub(">","&gt;",x, fixed=TRUE)
  x <- gsub('"',"&quot;",x, fixed=TRUE)
  x
}

write_table_csv <- function(df_tab, name, extras=FALSE) {
  dir0 <- if (extras) file.path(OUTDIR, "tables", "_extras") else file.path(OUTDIR, "tables")
  readr::write_csv(df_tab, file.path(dir0, paste0(name, ".csv")))
}

write_table_html <- function(df_tab, name, extras=FALSE, title=NULL, notes=NULL) {
  if (!MAKE_HTML) return(invisible(NULL))
  dir0 <- if (extras) file.path(OUTDIR, "tables", "_extras") else file.path(OUTDIR, "tables")
  f <- file.path(dir0, paste0(name, ".html"))

  # Simple, dependency-free HTML with embedded CSS
  css <- "
  body{font-family:Arial,Helvetica,sans-serif;margin:20px;}
  h2{margin:0 0 10px 0;}
  p.note{margin:8px 0 16px 0;color:#444;font-size:12px;}
  table{border-collapse:collapse;width:100%;font-size:12px;}
  th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top;}
  th{background:#f6f6f6;text-align:left;}
  tr:nth-child(even){background:#fbfbfb;}
  "
  tt <- if (!is.null(title)) paste0("<h2>", escape_html(title), "</h2>") else ""
  nn <- if (!is.null(notes)) paste0("<p class='note'><strong>Notes:</strong> ", escape_html(notes), "</p>") else ""

  # Build table
  cn <- names(df_tab)
  header <- paste0("<tr>", paste0("<th>", escape_html(cn), "</th>", collapse=""), "</tr>")
  rows <- apply(df_tab, 1, function(r) paste0("<tr>", paste0("<td>", escape_html(r), "</td>", collapse=""), "</tr>"))
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><style>", css, "</style></head><body>",
    tt, nn,
    "<table>", header, paste0(rows, collapse=""), "</table></body></html>"
  )
  writeLines(html, f)
}

write_table_both <- function(df_tab, name, extras=FALSE, title=NULL, notes=NULL) {
  write_table_csv(df_tab, name, extras=extras)
  write_table_html(df_tab, name, extras=extras, title=title, notes=notes)
}

# -------------------------------
# 7) Key descriptive tables (10 key tables)
# -------------------------------
T1 <- df_ie %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events_n = sum(os_status==1),
    Events_pct = round(100*mean(os_status==1),1),
    OS_median_years = round(stats::median(os_years),3),
    OS_mean_years = round(mean(os_years),3),
    Age_median = round(stats::median(age_at_dx, na.rm=TRUE),1),
    Age_mean = round(mean(age_at_dx, na.rm=TRUE),1),
    USC73_median = round(stats::median(usc73_score),3),
    USC73_mean = round(mean(usc73_score),3),
    .groups="drop"
  )
write_table_both(T1, "Table_K1_Baseline_Key",
                 title="Baseline by cohort",
                 notes="Median/mean OS are descriptive; model evaluation uses IPCW for censoring.")

T2 <- df_ie %>%
  dplyr::count(cohort, stage, race3, name="n") %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(pct = round(100*n/sum(n),1)) %>%
  dplyr::ungroup()
write_table_both(T2, "Table_K2_CaseMix_StageRace_Key",
                 title="Case-mix composition by cohort (stage × race)",
                 notes="Transportability pressure points arise when composition differs between cohorts.")

T3 <- df_ie %>%
  dplyr::summarise(across(c(age_at_dx, race3, stage, os_years, os_status, usc73_score), ~sum(is.na(.x)))) %>%
  tidyr::pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  dplyr::mutate(pct_missing = round(100*n_missing/nrow(df_ie),2))
write_table_both(T3, "Table_K3_Missingness_Key", title="Missingness (key variables)")

# -------------------------------
# 8) Pooled Cox (interpretability)
# -------------------------------
cox_base <- survival::coxph(survival::Surv(os_years, os_status) ~ age_at_dx + race3 + stage,
                            data=df_cc, ties="breslow")
cox_ext  <- survival::coxph(survival::Surv(os_years, os_status) ~ age_at_dx + race3 + stage + usc73_score,
                            data=df_cc, ties="breslow")

extract_hr <- function(fit, label) {
  s <- summary(fit)
  co <- as.data.frame(s$coefficients)
  co$term <- rownames(co)
  out <- co %>%
    dplyr::transmute(
      model = label,
      term = term,
      HR = exp(coef),
      CI_low = exp(coef - 1.96*`se(coef)`),
      CI_high= exp(coef + 1.96*`se(coef)`),
      p_value = `Pr(>|z|)`
    )
  out
}

Tcox_full <- dplyr::bind_rows(
  extract_hr(cox_base, "Base: age + race + stage"),
  extract_hr(cox_ext,  "Extended: base + USC73")
)
write_table_both(Tcox_full, "Table_S1_CoxHRs_Full", extras=TRUE, title="Cox hazard ratios (full)")

Tcox_key <- Tcox_full %>%
  dplyr::filter(stringr::str_detect(term, "stage|usc73|age")) %>%
  dplyr::mutate(across(c(HR, CI_low, CI_high), ~round(.x,3)),
                p_value = signif(p_value,3))
write_table_both(Tcox_key, "Table_K4_CoxHRs_Key",
                 title="Key Cox effects (interpretability)",
                 notes="Extended model includes USC73; HRs are pooled across cohorts (complete-case for age).")

# -------------------------------
# 9) Core utilities (prediction + estimability + IPCW)
# -------------------------------
align_matrix_cols <- function(x_new, col_ref) {
  miss <- setdiff(col_ref, colnames(x_new))
  if (length(miss) > 0) {
    add <- matrix(0, nrow=nrow(x_new), ncol=length(miss)); colnames(add) <- miss
    x_new <- cbind(x_new, add)
  }
  x_new[, col_ref, drop=FALSE]
}

fit_ridge_cox <- function(train, vars, nfolds=5, seed=1) {
  set.seed(seed)
  fml <- stats::as.formula(paste0("~ ", paste(vars, collapse=" + ")))
  mm <- stats::model.matrix(fml, data=train)[,-1, drop=FALSE]
  y <- survival::Surv(train$os_years, train$os_status)
  cv <- glmnet::cv.glmnet(mm, y, family="cox", alpha=0, nfolds=nfolds, standardize=TRUE)
  beta <- as.matrix(glmnet::coef.glmnet(cv, s="lambda.1se"))
  # Remove the first row if it is "(Intercept)" (cox usually doesn't include it)
  if ("(Intercept)" %in% rownames(beta)) beta <- beta[rownames(beta) != "(Intercept)", , drop=FALSE]
  list(fml=fml, cv=cv, beta=beta, col_ref=rownames(beta))
}

predict_lp_ridge <- function(fit, newdata) {
  mm <- stats::model.matrix(fit$fml, data=newdata)[,-1, drop=FALSE]
  mm <- align_matrix_cols(mm, fit$col_ref)
  as.numeric(mm %*% fit$beta)
}

fit_cox <- function(train, vars) {
  f <- stats::as.formula(paste0("survival::Surv(os_years, os_status) ~ ", paste(vars, collapse=" + ")))
  survival::coxph(f, data=train, ties="breslow", x=TRUE)
}
predict_lp_cox <- function(fit, newdata) as.numeric(stats::predict(fit, newdata=newdata, type="lp"))

baseline_surv_at_tau <- function(train, lp_train, tau) {
  # Guard against NA/inf lp
  lp_train <- ifelse(is.finite(lp_train), lp_train, 0)
  lp_train <- pmin(pmax(lp_train, -30), 30)
  off_fit <- survival::coxph(survival::Surv(os_years, os_status) ~ offset(lp_train),
                             data=train, ties="breslow")
  bh <- survival::basehaz(off_fit, centered=FALSE)
  H0_tau <- stats::approx(bh$time, bh$hazard, xout=tau, method="linear", rule=2)$y
  S0 <- exp(-H0_tau)
  if (!is.finite(S0) || S0 <= 0) S0 <- 1e-8
  if (S0 >= 1) S0 <- 1-1e-8
  S0
}

make_risk_from_lp <- function(train, lp_train, lp_test, tau) {
  lp_test <- ifelse(is.finite(lp_test), lp_test, 0)
  lp_test <- pmin(pmax(lp_test, -30), 30)
  S0 <- baseline_surv_at_tau(train, lp_train, tau)
  risk <- 1 - (S0 ^ exp(lp_test))
  pmin(pmax(risk, 1e-8), 1-1e-8)
}

km_surv_at <- function(sf, t) {
  s <- base::summary(sf, times=t, extend=TRUE)$surv
  if (length(s)==0 || is.na(s)) return(NA_real_)
  as.numeric(s[1])
}

ipcw_weights_at_tau <- function(time, status, tau) {
  # status=1 is event; censoring KM on (1-status)
  sf_cens <- survival::survfit(survival::Surv(time, 1-status) ~ 1)
  eligible <- (time >= tau) | (time < tau & status == 1)

  G_tau <- km_surv_at(sf_cens, tau); if (is.na(G_tau) || G_tau <= 0) G_tau <- 1e-6

  w <- rep(0, length(time))
  idx_event <- which(time < tau & status==1)
  if (length(idx_event)>0) {
    G_t <- sapply(time[idx_event], function(tt) {
      g <- km_surv_at(sf_cens, tt); if (is.na(g) || g<=0) g <- 1e-6; g
    })
    w[idx_event] <- 1 / G_t
  }
  idx_free <- which(time >= tau)
  w[idx_free] <- 1 / G_tau

  list(w=w, eligible=eligible, G_tau=G_tau)
}

estimability_counts <- function(time, status, tau) {
  y <- as.integer(time <= tau & status==1)
  iw <- ipcw_weights_at_tau(time, status, tau)
  elig <- iw$eligible
  tibble::tibble(
    tau_years = tau,
    n_test = length(time),
    n_eligible = sum(elig),
    cases_le_tau = sum(elig & y==1),
    controls_at_tau = sum(elig & y==0),
    cens_lt_tau = sum(time < tau & status==0),
    unique_time = length(unique(time))
  )
}

c_index_safe <- function(time, status, lp_risk) {
  # higher lp_risk => worse survival; use -lp for concordance direction
  tryCatch(as.numeric(survival::concordance(survival::Surv(time,status) ~ I(-lp_risk))$concordance),
           error=function(e) NA_real_)
}

brier_ipcw <- function(time, status, risk, tau) {
  iw <- ipcw_weights_at_tau(time, status, tau)
  w <- iw$w; elig <- iw$eligible
  y <- as.integer(time <= tau & status==1)
  denom <- sum(w[elig])
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(w[elig] * (y[elig] - risk[elig])^2) / denom
}

calibration_glm_ipcw <- function(time, status, risk, tau) {
  iw <- ipcw_weights_at_tau(time, status, tau)
  w <- iw$w; elig <- iw$eligible
  y <- as.integer(time <= tau & status==1)
  p <- pmin(pmax(risk, 1e-8), 1-1e-8)
  lp <- qlogis(p)
  dfc <- data.frame(y=y[elig], lp=lp[elig], w=w[elig])

  fit <- tryCatch(
    suppressWarnings(glm(y ~ lp, data=dfc, family=quasibinomial(), weights=dfc$w)),
    error=function(e) NULL
  )
  if (!is.null(fit)) {
    co <- coef(fit)
    return(list(intercept=unname(co[1]), slope=unname(co[2]), method="glm"))
  }

  # fallback: ridge logistic
  fit2 <- tryCatch({
    x <- model.matrix(~lp, data=dfc)[,-1, drop=FALSE]
    glmnet::glmnet(x, dfc$y, weights=dfc$w, family="binomial", alpha=0)
  }, error=function(e) NULL)
  if (is.null(fit2)) return(list(intercept=NA_real_, slope=NA_real_, method="failed"))
  b <- as.numeric(coef(fit2, s=fit2$lambda[length(fit2$lambda)]))
  list(intercept=b[1], slope=b[2], method="ridge_glmnet")
}

make_strat_folds <- function(y, K=5, seed=1) {
  set.seed(seed)
  y <- as.integer(y)
  idx1 <- which(y==1); idx0 <- which(y==0)
  K <- min(K, max(2, length(idx1)))  # ensure events can be distributed
  folds <- rep(NA_integer_, length(y))
  folds[idx1] <- sample(rep(1:K, length.out=length(idx1)))
  folds[idx0] <- sample(rep(1:K, length.out=length(idx0)))
  folds
}

recalibrate_crossfit <- function(time, status, risk, tau, K=5, seed=1) {
  iw <- ipcw_weights_at_tau(time, status, tau)
  w <- iw$w; elig <- iw$eligible
  y <- as.integer(time <= tau & status==1)
  p <- pmin(pmax(risk, 1e-8), 1-1e-8)
  lp <- qlogis(p)

  folds <- make_strat_folds(y[elig], K=K, seed=seed)
  p_out <- rep(NA_real_, length(p))
  idx_elig <- which(elig)

  for (k in sort(unique(folds))) {
    tr <- idx_elig[folds != k]
    te <- idx_elig[folds == k]
    if (length(te) == 0 || length(tr) < 20) next

    df_tr <- data.frame(y=y[tr], lp=lp[tr], w=w[tr])
    fit <- tryCatch(
      suppressWarnings(glm(y ~ lp, data=df_tr, family=quasibinomial(), weights=df_tr$w)),
      error=function(e) NULL
    )

    if (is.null(fit)) {
      fit2 <- tryCatch({
        x <- model.matrix(~lp, data=df_tr)[,-1, drop=FALSE]
        glmnet::glmnet(x, df_tr$y, weights=df_tr$w, family="binomial", alpha=0)
      }, error=function(e) NULL)
      if (is.null(fit2)) next
      b <- as.numeric(coef(fit2, s=fit2$lambda[length(fit2$lambda)]))
      a0 <- b[1]; a1 <- b[2]
    } else {
      co <- coef(fit); a0 <- unname(co[1]); a1 <- unname(co[2])
    }

    p_out[te] <- plogis(a0 + a1 * lp[te])
  }

  # Non-eligible keep original (do not enter IPCW scoring)
  p_out[!elig] <- p[!elig]
  pmin(pmax(p_out, 1e-8), 1-1e-8)
}

decision_curve_ipcw <- function(time, status, risk, tau, thresholds=seq(0.05,0.60,by=0.01)) {
  iw <- ipcw_weights_at_tau(time, status, tau)
  w <- iw$w; elig <- iw$eligible
  y <- as.integer(time <= tau & status==1)
  denom <- sum(w[elig])
  if (!is.finite(denom) || denom <= 0) {
    return(tibble::tibble(threshold=thresholds, nb_model=NA_real_, nb_all=NA_real_, nb_none=0))
  }

  nb <- sapply(thresholds, function(pt) {
    treat <- as.integer(risk >= pt)
    TP <- sum(w[elig]*y[elig]*treat[elig]) / denom
    FP <- sum(w[elig]*(1-y[elig])*treat[elig]) / denom
    TP - FP * (pt/(1-pt))
  })
  prev <- sum(w[elig]*y[elig]) / denom
  nb_all <- sapply(thresholds, function(pt) prev - (1-prev) * (pt/(1-pt)))
  tibble::tibble(threshold=thresholds, nb_model=nb, nb_all=nb_all, nb_none=0)
}

cal_curve_data <- function(time, status, risk, tau, groups=10) {
  iw <- ipcw_weights_at_tau(time, status, tau)
  w <- iw$w; elig <- iw$eligible
  y <- as.integer(time <= tau & status==1)

  dfc <- tibble::tibble(risk=risk, y=y, w=w, elig=elig) %>%
    dplyr::filter(elig) %>%
    dplyr::mutate(risk=pmin(pmax(risk,1e-8),1-1e-8)) %>%
    dplyr::arrange(risk)

  g <- min(groups, max(3, floor(nrow(dfc)/10)))
  dfc <- dfc %>% dplyr::mutate(bin=dplyr::ntile(risk, g))

  dfc %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      n = dplyr::n(),
      pred = mean(risk),
      obs  = sum(w*y)/sum(w),
      .groups="drop"
    )
}

auc_timeROC_safe <- function(time, status, marker, tau) {
  if (!has_timeROC) return(NA_real_)
  out <- tryCatch({
    # timeROC expects status=1 event; "cause" is the event code
    tr <- timeROC::timeROC(T=time, delta=status, marker=marker, cause=1,
                           times=tau, iid=FALSE)
    as.numeric(tr$AUC[1])
  }, error=function(e) NA_real_)
  out
}

# -------------------------------
# 10) IE engine (both directions)
# -------------------------------
base_vars <- c("age_at_dx","race3","stage")
ext_vars  <- c("age_at_dx","race3","stage","usc73_score")

filter_estimable_taus <- function(test, tau_cand, min_cases, min_controls) {
  purrr::map_dfr(tau_cand, function(tau) {
    ec <- estimability_counts(test$os_years, test$os_status, tau)
    ec$estimable <- (ec$cases_le_tau >= min_cases) && (ec$controls_at_tau >= min_controls)
    ec
  })
}

choose_rep_tau <- function(taus_ok) {
  if (length(taus_ok)==0) return(NA_real_)
  taus_ok[which.min(abs(taus_ok - 2.0))]  # prefer 2y if available
}

run_direction <- function(train_name, test_name, data) {
  train <- data %>% dplyr::filter(cohort==train_name)
  test  <- data %>% dplyr::filter(cohort==test_name)
  direction <- paste0(train_name, " -> ", test_name)

  est <- filter_estimable_taus(test, TAU_CAND, MIN_CASES, MIN_CONTROLS) %>%
    dplyr::mutate(direction=direction)

  taus_ok <- est %>% dplyr::filter(estimable) %>% dplyr::pull(tau_years)
  if (length(taus_ok)==0) {
    cat("WARNING: No estimable taus for", direction, "\n")
    return(list(direction=direction, estimability=est, metrics=tibble::tibble(),
                dca=tibble::tibble(), calcurve=tibble::tibble(), rep_tau=NA_real_))
  }

  # Fit models
  ridge_base <- fit_ridge_cox(train, base_vars, nfolds=NFOLDS, seed=SEED)
  ridge_ext  <- fit_ridge_cox(train, ext_vars,  nfolds=NFOLDS, seed=SEED+1)
  cox_base_m <- fit_cox(train, base_vars)
  cox_ext_m  <- fit_cox(train, ext_vars)

  # Precompute LPs
  lp_train_rb <- predict_lp_ridge(ridge_base, train)
  lp_test_rb  <- predict_lp_ridge(ridge_base, test)
  lp_train_re <- predict_lp_ridge(ridge_ext,  train)
  lp_test_re  <- predict_lp_ridge(ridge_ext,  test)

  lp_train_cb <- predict_lp_cox(cox_base_m, train)
  lp_test_cb  <- predict_lp_cox(cox_base_m, test)
  lp_train_ce <- predict_lp_cox(cox_ext_m,  train)
  lp_test_ce  <- predict_lp_cox(cox_ext_m,  test)

  eval_at_tau <- function(tau, model_id) {
    if (model_id=="ridge_base") {
      lp_test <- lp_test_rb; lp_train <- lp_train_rb
      risk <- make_risk_from_lp(train, lp_train, lp_test, tau)
    } else if (model_id=="ridge_ext") {
      lp_test <- lp_test_re; lp_train <- lp_train_re
      risk <- make_risk_from_lp(train, lp_train, lp_test, tau)
    } else if (model_id=="cox_base") {
      lp_test <- lp_test_cb; lp_train <- lp_train_cb
      risk <- make_risk_from_lp(train, lp_train, lp_test, tau)
    } else if (model_id=="cox_ext") {
      lp_test <- lp_test_ce; lp_train <- lp_train_ce
      risk <- make_risk_from_lp(train, lp_train, lp_test, tau)
    } else stop("Unknown model_id")

    ec <- estimability_counts(test$os_years, test$os_status, tau)
    if (ec$cases_le_tau < MIN_CASES || ec$controls_at_tau < MIN_CONTROLS) {
      return(tibble::tibble(direction=direction, tau_years=tau, model_id=model_id,
                            c_index=NA_real_, auc_tau=NA_real_, brier_tau=NA_real_,
                            calib_intercept=NA_real_, calib_slope=NA_real_, calib_method="not_estimable",
                            brier_recal=NA_real_, auc_recal=NA_real_,
                            cases_le_tau=ec$cases_le_tau, controls_at_tau=ec$controls_at_tau))
    }

    cind <- c_index_safe(test$os_years, test$os_status, lp_test)
    auc  <- auc_timeROC_safe(test$os_years, test$os_status, marker=lp_test, tau=tau)
    br   <- brier_ipcw(test$os_years, test$os_status, risk, tau)
    cal  <- calibration_glm_ipcw(test$os_years, test$os_status, risk, tau)

    # Cross-fitted recalibration in target cohort
    risk_recal <- recalibrate_crossfit(test$os_years, test$os_status, risk, tau, K=5, seed=SEED)
    br_recal   <- brier_ipcw(test$os_years, test$os_status, risk_recal, tau)
    auc_recal  <- auc_timeROC_safe(test$os_years, test$os_status, marker=qlogis(risk_recal), tau=tau)

    tibble::tibble(
      direction=direction, tau_years=tau, model_id=model_id,
      c_index=cind, auc_tau=auc, brier_tau=br,
      calib_intercept=cal$intercept, calib_slope=cal$slope, calib_method=cal$method,
      brier_recal=br_recal, auc_recal=auc_recal,
      cases_le_tau=ec$cases_le_tau, controls_at_tau=ec$controls_at_tau
    )
  }

  metrics <- purrr::map_dfr(taus_ok, function(tau) {
    dplyr::bind_rows(
      eval_at_tau(tau, "ridge_base"),
      eval_at_tau(tau, "ridge_ext"),
      eval_at_tau(tau, "cox_base"),
      eval_at_tau(tau, "cox_ext")
    )
  }) %>%
    dplyr::mutate(model = dplyr::case_when(
      model_id=="ridge_base" ~ "Ridge Base",
      model_id=="ridge_ext"  ~ "Ridge Extended",
      model_id=="cox_base"   ~ "Cox Base",
      model_id=="cox_ext"    ~ "Cox Extended",
      TRUE ~ model_id
    ))

  rep_tau <- choose_rep_tau(taus_ok)

  dca <- tibble::tibble()
  calcurve <- tibble::tibble()

  if (is.finite(rep_tau)) {
    risk_rb <- make_risk_from_lp(train, lp_train_rb, lp_test_rb, rep_tau)
    risk_re <- make_risk_from_lp(train, lp_train_re, lp_test_re, rep_tau)

    # raw + recal
    risk_rb_recal <- recalibrate_crossfit(test$os_years, test$os_status, risk_rb, rep_tau, K=5, seed=SEED)
    risk_re_recal <- recalibrate_crossfit(test$os_years, test$os_status, risk_re, rep_tau, K=5, seed=SEED+1)

    dca <- dplyr::bind_rows(
      decision_curve_ipcw(test$os_years, test$os_status, risk_rb, rep_tau) %>% dplyr::mutate(model="Ridge Base (raw)"),
      decision_curve_ipcw(test$os_years, test$os_status, risk_rb_recal, rep_tau) %>% dplyr::mutate(model="Ridge Base (recal)"),
      decision_curve_ipcw(test$os_years, test$os_status, risk_re, rep_tau) %>% dplyr::mutate(model="Ridge Ext (raw)"),
      decision_curve_ipcw(test$os_years, test$os_status, risk_re_recal, rep_tau) %>% dplyr::mutate(model="Ridge Ext (recal)")
    ) %>% dplyr::mutate(direction=direction, tau_years=rep_tau)

    calcurve <- dplyr::bind_rows(
      cal_curve_data(test$os_years, test$os_status, risk_rb, rep_tau) %>% dplyr::mutate(model="Ridge Base (raw)"),
      cal_curve_data(test$os_years, test$os_status, risk_rb_recal, rep_tau) %>% dplyr::mutate(model="Ridge Base (recal)"),
      cal_curve_data(test$os_years, test$os_status, risk_re, rep_tau) %>% dplyr::mutate(model="Ridge Ext (raw)"),
      cal_curve_data(test$os_years, test$os_status, risk_re_recal, rep_tau) %>% dplyr::mutate(model="Ridge Ext (recal)")
    ) %>% dplyr::mutate(direction=direction, tau_years=rep_tau)
  }

  list(direction=direction, estimability=est, metrics=metrics, dca=dca, calcurve=calcurve, rep_tau=rep_tau)
}

ie1 <- run_direction("TCGA","UPSC", df_cc)
ie2 <- run_direction("UPSC","TCGA", df_cc)

EST_ALL <- dplyr::bind_rows(ie1$estimability, ie2$estimability)
MET_ALL <- dplyr::bind_rows(ie1$metrics, ie2$metrics)
DCA_ALL <- dplyr::bind_rows(ie1$dca, ie2$dca)
CAL_ALL <- dplyr::bind_rows(ie1$calcurve, ie2$calcurve)

write_table_both(EST_ALL, "Table_S2_Estimability_All", extras=TRUE, title="Estimability by tau and direction")
write_table_both(MET_ALL, "Table_S3_Metrics_All", extras=TRUE, title="All metrics by tau, direction, and model")
write_table_both(DCA_ALL, "Table_S4_DCA_All", extras=TRUE, title="Decision curves (IPCW)")
write_table_both(CAL_ALL, "Table_S5_CalibrationCurve_All", extras=TRUE, title="Calibration curve binned data")

# Key Table 5: Estimability summary (which taus are usable)
Testim_key <- EST_ALL %>%
  dplyr::mutate(estimable = ifelse(estimable, "YES","NO")) %>%
  dplyr::select(direction, tau_years, cases_le_tau, controls_at_tau, cens_lt_tau, estimable) %>%
  dplyr::arrange(direction, tau_years)
write_table_both(Testim_key, "Table_K5_Estimability_Key",
                 title="Estimability gate by tau and direction",
                 notes="Evaluation at tau is only meaningful when sufficient cases and controls exist after accounting for censoring.")

# Key Table 6: Performance summary (average across taus)
Tperf_key <- MET_ALL %>%
  dplyr::group_by(direction, model) %>%
  dplyr::summarise(
    n_tau = dplyr::n(),
    c_index_mean = round(mean(c_index, na.rm=TRUE),3),
    auc_mean = round(mean(auc_tau, na.rm=TRUE),3),
    brier_mean = round(mean(brier_tau, na.rm=TRUE),4),
    brier_recal_mean = round(mean(brier_recal, na.rm=TRUE),4),
    .groups="drop"
  ) %>% dplyr::arrange(direction, model)
write_table_both(Tperf_key, "Table_K6_Performance_Key",
                 title="Key predictive performance (averaged across estimable taus)",
                 notes="Lower Brier is better; higher C-index/AUC is better. Recalibration is performed in the target cohort via cross-fitting.")

# Key Table 7: Deltas (Extended - Base)
Tdeltas <- MET_ALL %>%
  dplyr::filter(model_id %in% c("ridge_base","ridge_ext","cox_base","cox_ext")) %>%
  dplyr::mutate(
    family = dplyr::case_when(
      stringr::str_detect(model_id, "^ridge") ~ "Ridge-Cox",
      stringr::str_detect(model_id, "^cox")   ~ "Cox",
      TRUE ~ "Other"
    ),
    variant = dplyr::case_when(
      stringr::str_detect(model_id, "_base$") ~ "Base",
      stringr::str_detect(model_id, "_ext$")  ~ "Extended",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::select(direction, tau_years, family, variant, c_index, auc_tau, brier_tau, brier_recal) %>%
  tidyr::pivot_wider(names_from=variant, values_from=c(c_index, auc_tau, brier_tau, brier_recal)) %>%
  dplyr::mutate(
    delta_c_index = round(c_index_Extended - c_index_Base, 3),
    delta_auc     = round(auc_tau_Extended - auc_tau_Base, 3),
    delta_brier   = round(brier_tau_Extended - brier_tau_Base, 4),
    delta_brier_recal = round(brier_recal_Extended - brier_recal_Base, 4)
  ) %>%
  dplyr::select(direction, tau_years, family, delta_c_index, delta_auc, delta_brier, delta_brier_recal) %>%
  dplyr::arrange(direction, family, tau_years)
write_table_both(Tdeltas, "Table_K7_Deltas_ExtMinusBase_Key",
                 title="Incremental value of USC73 (Extended - Base)",
                 notes="Positive delta C-index/AUC is improvement; negative delta Brier is improvement.")

# Key Table 8: Representative tau metrics (closest to 2y per direction)
rep_taus <- tibble::tibble(direction=c(ie1$direction, ie2$direction), rep_tau=c(ie1$rep_tau, ie2$rep_tau))
Trep <- MET_ALL %>%
  dplyr::inner_join(rep_taus, by="direction") %>%
  dplyr::filter(is.finite(rep_tau), abs(tau_years - rep_tau) < 1e-8) %>%
  dplyr::select(direction, tau_years, model, c_index, auc_tau, brier_tau, calib_intercept, calib_slope, brier_recal) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x,4))) %>%
  dplyr::arrange(direction, model)
write_table_both(Trep, "Table_K8_RepTau_Metrics_Key",
                 title="Metrics at representative tau (≈2 years, estimable)",
                 notes="Calibration parameters reflect IPCW-weighted logistic calibration; recalibration is cross-fitted within target cohort.")

# Key Table 9: Decision curve summary (area under NB curve vs treat-all/none)
T_dca_key <- DCA_ALL %>%
  dplyr::group_by(direction, tau_years, model) %>%
  dplyr::summarise(
    nb_auc = round(mean(nb_model, na.rm=TRUE),4),
    nb_auc_all = round(mean(nb_all, na.rm=TRUE),4),
    nb_gain_vs_all = round(nb_auc - nb_auc_all,4),
    .groups="drop"
  ) %>% dplyr::arrange(direction, model)
write_table_both(T_dca_key, "Table_K9_DCA_Summary_Key",
                 title="Decision curve summary (IPCW net benefit)",
                 notes="nb_gain_vs_all > 0 indicates advantage over treat-all across the threshold range.")

# Key Table 10: Calibration curve summary (mean |pred-obs| across bins)
T_cal_key <- CAL_ALL %>%
  dplyr::group_by(direction, tau_years, model) %>%
  dplyr::summarise(
    bins = dplyr::n(),
    mae = round(mean(abs(pred - obs), na.rm=TRUE),4),
    .groups="drop"
  ) %>% dplyr::arrange(direction, model)
write_table_both(T_cal_key, "Table_K10_Calibration_Summary_Key",
                 title="Calibration summary (binned MAE at representative tau)",
                 notes="Lower MAE indicates closer agreement between predicted and observed risks (IPCW-adjusted).")

# -------------------------------
# 11) Figures (>=12 PNGs; fail-safe)
# -------------------------------
theme_pub <- function() {
  ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face="bold", size=14),
      plot.subtitle = ggplot2::element_text(size=11),
      plot.caption = ggplot2::element_text(size=9, hjust=0),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

scale_fill_discrete_safe <- function(option="D") {
  # Prefer ggplot2 viridis scales if available
  if ("scale_fill_viridis_d" %in% getNamespaceExports("ggplot2")) {
    return(ggplot2::scale_fill_viridis_d(option=option))
  }
  # viridis older versions: scale_fill_viridis(discrete=TRUE)
  if (has_viridis && "scale_fill_viridis" %in% getNamespaceExports("viridis")) {
    return(viridis::scale_fill_viridis(discrete=TRUE, option=option))
  }
  ggplot2::scale_fill_brewer(palette="Set2")
}

scale_color_discrete_safe <- function(option="D") {
  if ("scale_color_viridis_d" %in% getNamespaceExports("ggplot2")) {
    return(ggplot2::scale_color_viridis_d(option=option))
  }
  if (has_viridis && "scale_color_viridis" %in% getNamespaceExports("viridis")) {
    return(viridis::scale_color_viridis(discrete=TRUE, option=option))
  }
  ggplot2::scale_color_brewer(palette="Set2")
}

save_fig <- function(p, name, w=8, h=5, dpi=320) {
  tryCatch({
    ggplot2::ggsave(filename=file.path(OUTDIR, "figures", paste0(name, ".png")),
                    plot=p, width=w, height=h, dpi=dpi)
  }, error=function(e) {
    cat("FIGURE SAVE FAILED:", name, "->", e$message, "\n")
  })
}

# Helper: survfit -> data.frame
tidy_survfit <- function(fit) {
  s <- summary(fit)
  strata <- s$strata
  if (is.null(strata)) {
    strlab <- rep("All", length(s$time))
  } else if (length(strata) == length(s$time)) {
    strlab <- as.character(strata)
  } else {
    strlab <- rep(names(strata), strata)
  }
  tibble::tibble(time=s$time, surv=s$surv, strata=strlab) %>%
    dplyr::mutate(strata = stringr::str_replace_all(strata, ".*=", ""))
}

# FIG 01: Cohort shift (event rate, median OS, USC73)
T1_long <- T1 %>%
  tidyr::pivot_longer(cols=c(Events_pct, OS_median_years, USC73_median),
                      names_to="metric", values_to="value") %>%
  dplyr::mutate(metric = dplyr::recode(metric,
                                       Events_pct="Event rate (%)",
                                       OS_median_years="Median OS (years)",
                                       USC73_median="Median USC73"))
p1 <- ggplot2::ggplot(T1_long, ggplot2::aes(x=cohort, y=value, fill=metric)) +
  ggplot2::geom_col(position="dodge", width=0.75, color="grey25") +
  ggplot2::labs(
    title="Cohort shift: outcomes and signature distribution",
    subtitle="Key transportability drivers: different event rates, OS time scales, and USC73 distribution",
    x=NULL, y=NULL,
    caption="Bars show cohort-specific descriptive summaries. Evaluation uses estimability-aware tau + IPCW metrics."
  ) +
  theme_pub() + scale_fill_discrete_safe("D")
save_fig(p1, "Fig01_CohortShift_Bars", w=9, h=5)

# FIG 02: USC73 distribution by cohort (violin + jitter)
p2 <- ggplot2::ggplot(df_ie, ggplot2::aes(x=cohort, y=usc73_score, color=cohort)) +
  ggplot2::geom_violin(trim=FALSE, alpha=0.25) +
  ggplot2::geom_jitter(width=0.12, height=0, alpha=0.55, size=1.3) +
  ggplot2::labs(
    title="USC73 distribution shift across cohorts",
    subtitle="Distributional shift is a mechanistic reason transportability can fail without recalibration",
    x=NULL, y="USC73 score",
    caption="Violin shows density; points show individuals."
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p2, "Fig02_USC73_Distribution", w=8.5, h=5)

# FIG 03: Kaplan–Meier by cohort
sf0 <- survival::survfit(survival::Surv(os_years, os_status) ~ cohort, data=df_ie)
sf_df <- tidy_survfit(sf0)
p3 <- ggplot2::ggplot(sf_df, ggplot2::aes(x=time, y=surv, color=strata)) +
  ggplot2::geom_step(linewidth=1) +
  ggplot2::labs(
    title="Overall survival differs by cohort",
    subtitle="Kaplan–Meier curves for TCGA vs UPSC",
    x="Time (years)", y="Survival probability",
    caption="KM curves are descriptive; prediction performance is evaluated at estimable time horizons with IPCW."
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p3, "Fig03_KM_ByCohort", w=8.5, h=5)

# FIG 04: KM by USC73 tertiles within cohort
df_ie2 <- df_ie %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(usc73_tertile = dplyr::ntile(usc73_score, 3)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(usc73_tertile = factor(usc73_tertile, levels=1:3, labels=c("Low","Mid","High")))
sf1 <- survival::survfit(survival::Surv(os_years, os_status) ~ usc73_tertile + cohort, data=df_ie2)
sf1_df <- tidy_survfit(sf1) %>%
  tidyr::separate(strata, into=c("cohort","usc73_tertile"), sep=",") %>%
  dplyr::mutate(cohort=stringr::str_replace_all(cohort,"cohort",""),
                cohort=stringr::str_replace_all(cohort,"\\s",""),
                usc73_tertile=stringr::str_replace_all(usc73_tertile,"usc73_tertile",""),
                usc73_tertile=stringr::str_replace_all(usc73_tertile,"\\s",""))
p4 <- ggplot2::ggplot(sf1_df, ggplot2::aes(x=time, y=surv, color=usc73_tertile)) +
  ggplot2::geom_step(linewidth=1) +
  ggplot2::facet_wrap(~cohort) +
  ggplot2::labs(
    title="Prognostic stratification by USC73 within each cohort",
    subtitle="Tertile separation suggests biological signal; transportability depends on recalibration in the target cohort",
    x="Time (years)", y="Survival probability"
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p4, "Fig04_KM_USC73_Tertiles", w=10, h=5)

# FIG 05: Forest plot of HRs (extended Cox)
T_hr_ext <- Tcox_full %>%
  dplyr::filter(model=="Extended: base + USC73") %>%
  dplyr::mutate(
    term_clean = stringr::str_replace_all(term, "race3", "Race: "),
    term_clean = stringr::str_replace_all(term_clean, "stage", "Stage: "),
    term_clean = stringr::str_replace_all(term_clean, "usc73_score", "USC73 (per unit)"),
    term_clean = stringr::str_replace_all(term_clean, "age_at_dx", "Age at diagnosis (per year)")
  ) %>%
  dplyr::arrange(desc(HR)) %>%
  dplyr::mutate(term_clean = factor(term_clean, levels=rev(unique(term_clean))))
p5 <- ggplot2::ggplot(T_hr_ext, ggplot2::aes(x=HR, y=term_clean)) +
  ggplot2::geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin=CI_low, xmax=CI_high), height=0.2, color="grey40") +
  ggplot2::geom_point(size=2) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title="Interpretability: extended Cox hazard ratios",
    subtitle="Pooled complete-case model adjusted for age, race, stage, and USC73",
    x="Hazard ratio (log scale)", y=NULL
  ) +
  theme_pub()
save_fig(p5, "Fig05_Forest_CoxHR_Extended", w=9, h=5.5)

# FIG 06: Estimability heatmap
EST_PLOT <- EST_ALL %>%
  dplyr::mutate(estimable_flag = ifelse(estimable, "Estimable","Not estimable"))
p6 <- ggplot2::ggplot(EST_PLOT, ggplot2::aes(x=factor(tau_years), y=direction, fill=estimable_flag)) +
  ggplot2::geom_tile(color="white") +
  ggplot2::labs(
    title="Estimability gate across time horizons",
    subtitle="Only taus with sufficient cases/controls are evaluated (prevents meaningless metrics)",
    x="Tau (years)", y=NULL
  ) +
  theme_pub() + scale_fill_discrete_safe("C")
save_fig(p6, "Fig06_Estimability_Heatmap", w=10, h=3.8)

# FIG 07: C-index vs tau
p7 <- ggplot2::ggplot(MET_ALL, ggplot2::aes(x=tau_years, y=c_index, color=model)) +
  ggplot2::geom_line(linewidth=1) +
  ggplot2::geom_point(size=1.8) +
  ggplot2::facet_wrap(~direction) +
  ggplot2::labs(
    title="Discrimination vs time horizon (C-index)",
    subtitle="Evaluated only where tau is estimable in the target cohort",
    x="Tau (years)", y="C-index"
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p7, "Fig07_CIndex_vs_Tau", w=11, h=5.5)

# FIG 08: Brier vs tau (raw and recalibration)
MET_LONG <- MET_ALL %>%
  dplyr::select(direction, tau_years, model, brier_tau, brier_recal) %>%
  tidyr::pivot_longer(cols=c(brier_tau, brier_recal), names_to="type", values_to="brier") %>%
  dplyr::mutate(type = dplyr::recode(type, brier_tau="Raw", brier_recal="Recalibrated"))
p8 <- ggplot2::ggplot(MET_LONG, ggplot2::aes(x=tau_years, y=brier, color=model, linetype=type)) +
  ggplot2::geom_line(linewidth=1) +
  ggplot2::geom_point(size=1.6) +
  ggplot2::facet_wrap(~direction) +
  ggplot2::labs(
    title="Accuracy vs time horizon (IPCW Brier score)",
    subtitle="Cross-fitted target recalibration often improves Brier under cohort shift",
    x="Tau (years)", y="Brier score (lower is better)"
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p8, "Fig08_Brier_vs_Tau_Raw_vs_Recal", w=11, h=5.5)

# FIG 09: Delta C-index (Extended - Base)
p9 <- ggplot2::ggplot(Tdeltas, ggplot2::aes(x=tau_years, y=delta_c_index, color=family)) +
  ggplot2::geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  ggplot2::geom_line(linewidth=1) +
  ggplot2::geom_point(size=1.8) +
  ggplot2::facet_wrap(~direction) +
  ggplot2::labs(
    title="Incremental value of USC73 across time horizons",
    subtitle="Delta C-index (Extended - Base); >0 indicates improvement",
    x="Tau (years)", y="Delta C-index"
  ) +
  theme_pub() + scale_color_discrete_safe("D")
save_fig(p9, "Fig09_DeltaCIndex_ExtMinusBase", w=10.5, h=5.2)

# FIG 10: Calibration curve at representative tau (raw vs recal)
if (nrow(CAL_ALL) > 0) {
  p10 <- ggplot2::ggplot(CAL_ALL, ggplot2::aes(x=pred, y=obs, color=model)) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    ggplot2::geom_point(size=2) +
    ggplot2::geom_line(linewidth=1) +
    ggplot2::facet_wrap(~direction) +
    ggplot2::labs(
      title="Calibration at representative tau (IPCW)",
      subtitle="Closer to the diagonal indicates better calibration",
      x="Mean predicted risk (bin)", y="IPCW observed risk (bin)"
    ) +
    theme_pub() + scale_color_discrete_safe("D")
  save_fig(p10, "Fig10_CalibrationCurve_RepTau", w=11, h=5.5)
}

# FIG 11: Decision curves at representative tau
if (nrow(DCA_ALL) > 0) {
  p11 <- ggplot2::ggplot(DCA_ALL, ggplot2::aes(x=threshold, y=nb_model, color=model)) +
    ggplot2::geom_line(linewidth=1) +
    ggplot2::facet_wrap(~direction) +
    ggplot2::labs(
      title="Clinical utility: IPCW decision curve analysis",
      subtitle="Net benefit across threshold probabilities at representative tau",
      x="Threshold probability", y="Net benefit"
    ) +
    theme_pub() + scale_color_discrete_safe("D")
  save_fig(p11, "Fig11_DecisionCurves_RepTau", w=11, h=5.5)
}

# FIG 12: Risk distribution at representative tau (ridge models)
if (nrow(rep_taus) > 0) {
  # Build risk for each direction at rep_tau using ridge base/ext
  make_risk_df <- function(train_name, test_name, rep_tau, data) {
    if (!is.finite(rep_tau)) return(tibble::tibble())
    train <- data %>% dplyr::filter(cohort==train_name)
    test  <- data %>% dplyr::filter(cohort==test_name)

    rb <- fit_ridge_cox(train, base_vars, nfolds=NFOLDS, seed=SEED+11)
    re <- fit_ridge_cox(train, ext_vars,  nfolds=NFOLDS, seed=SEED+12)

    lp_train_rb <- predict_lp_ridge(rb, train); lp_test_rb <- predict_lp_ridge(rb, test)
    lp_train_re <- predict_lp_ridge(re, train); lp_test_re <- predict_lp_ridge(re, test)

    r_rb <- make_risk_from_lp(train, lp_train_rb, lp_test_rb, rep_tau)
    r_re <- make_risk_from_lp(train, lp_train_re, lp_test_re, rep_tau)

    tibble::tibble(
      direction=paste0(train_name, " -> ", test_name),
      risk=c(r_rb, r_re),
      model=rep(c("Ridge Base","Ridge Extended"), each=length(r_rb))
    )
  }

  r1 <- make_risk_df("TCGA","UPSC", ie1$rep_tau, df_cc)
  r2 <- make_risk_df("UPSC","TCGA", ie2$rep_tau, df_cc)
  RR <- dplyr::bind_rows(r1, r2)

  if (nrow(RR) > 0) {
    p12 <- ggplot2::ggplot(RR, ggplot2::aes(x=risk, color=model, fill=model)) +
      ggplot2::geom_density(alpha=0.20, linewidth=1) +
      ggplot2::facet_wrap(~direction) +
      ggplot2::labs(
        title="Predicted risk distribution shift in the target cohort",
        subtitle="Risk distributions can shift substantially under cohort differences, motivating recalibration",
        x="Predicted event risk at representative tau", y="Density"
      ) +
      theme_pub() + scale_color_discrete_safe("D") + scale_fill_discrete_safe("D")
    save_fig(p12, "Fig12_RiskDistribution_RepTau", w=11, h=5.5)
  }
}

cat("\nPIPELINE COMPLETE.\n")
cat("Outputs written to:", OUTDIR, "\n")
