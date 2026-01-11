#!/usr/bin/env Rscript
# ============================================================
# 01_profile_usc73_excel.R
#
# Excel-first profiling + QC for:
#   Supplementary_Expression_phenotype_data.xlsx
#
# What it does:
#   - Enumerates Excel sheets, reads each, and scores likelihood of being the
#     merged "expression + phenotype" table (clinical keywords + many columns)
#   - Profiles the selected sheet (or a user-specified sheet)
#   - Outputs data dictionary, missingness, summaries, and QC flags
#
# Run:
#   Rscript 01_profile_usc73_excel.R --input "Supplementary_Expression_phenotype_data.xlsx" --outdir "qc_usc73"
#
# Optional:
#   --sheet "SheetName"     # force a sheet
#   --plots true|false      # default true
#   --guess_max 5000        # readxl type-guessing depth
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c("readxl", "readr", "dplyr", "tibble", "stringr", "janitor", "ggplot2", "tidyr", "purrr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  lapply(pkgs, library, character.only = TRUE)
})

# -------------------------------
# Minimal CLI parsing
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}
INPUT    <- get_arg("--input", default = "Supplementary_Expression_phenotype_data.xlsx")
OUTDIR   <- get_arg("--outdir", default = "qc_output")
SHEET_IN <- get_arg("--sheet", default = NA)
PLOTS    <- tolower(get_arg("--plots", default = "true")) %in% c("true","t","1","yes","y")
GUESSMAX <- suppressWarnings(as.integer(get_arg("--guess_max", default = "5000")))
if (is.na(GUESSMAX) || GUESSMAX < 1000) GUESSMAX <- 5000

if (!file.exists(INPUT)) stop("Input Excel file not found: ", INPUT, call. = FALSE)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (PLOTS) dir.create(file.path(OUTDIR, "plots"), showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Helpers
# -------------------------------
clean_names_unique <- function(x) make.unique(janitor::make_clean_names(x), sep = "__dup__")

# Clinical keyword heuristics (after clean_names)
clinical_regex <- "(^age$|age_at_dx|race|ethnic|stage|grade|os_years|os_status|overall_survival|survival|status|event|death|treat|therapy|response|usc73)"

score_sheet <- function(df) {
  cn <- janitor::make_clean_names(names(df))
  clin_hits <- sum(stringr::str_detect(cn, clinical_regex))
  # Prefer wide tables (genes) but require clinical hits
  # Weight clinical hits heavily (so phenotype-only or gene-only sheets don't win accidentally)
  score <- clin_hits * 1000 + ncol(df)
  tibble::tibble(clin_hits = clin_hits, ncol = ncol(df), nrow = nrow(df), score = score)
}

infer_role <- function(varname) {
  v <- tolower(varname)
  if (v %in% c("age_at_dx", "age", "race", "os_years", "os_status", "stage", "usc73")) return("clinical_core")
  if (grepl("^os_|survival|time|status|event|death", v)) return("clinical_outcome")
  if (grepl("age|race|ethnic|stage|grade|treat|therapy|response", v)) return("clinical_covariate")
  return("molecular_feature")
}

is_date_like <- function(x) {
  if (!is.character(x)) return(FALSE)
  xx <- x[!is.na(x)]
  if (length(xx) == 0) return(FALSE)
  samp <- xx[1:min(50, length(xx))]
  any(grepl("^\\d{4}-\\d{2}-\\d{2}$", samp)) || any(grepl("^\\d{2}/\\d{2}/\\d{4}$", samp))
}

infer_type <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) return("date_like")
  
  if (is.character(x)) {
    if (is_date_like(x)) return("date_like")
    xn <- suppressWarnings(as.numeric(x))
    ok <- sum(!is.na(xn)) / max(1, sum(!is.na(x)))
    if (ok >= 0.9 && sum(!is.na(x)) > 0) {
      if (all(abs(xn[!is.na(xn)] - round(xn[!is.na(xn)])) < 1e-8)) return("integer_like")
      return("numeric")
    }
    nuniq <- length(unique(x[!is.na(x)]))
    if (nuniq <= 20) return("categorical")
    return("text")
  }
  
  if (is.logical(x)) return("binary")
  if (is.integer(x)) {
    vals <- unique(x[!is.na(x)])
    if (length(vals) <= 2) return("binary")
    return("integer_like")
  }
  if (is.numeric(x)) {
    vals <- unique(x[!is.na(x)])
    if (length(vals) <= 2) return("binary")
    if (all(abs(x[!is.na(x)] - round(x[!is.na(x)])) < 1e-8)) return("integer_like")
    return("numeric")
  }
  return("text")
}

# Identical-column detection (robust across types)
col_identical <- function(a, b) identical(as.character(a), as.character(b))

profile_df <- function(df, outdir, plots = TRUE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (plots) dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Preserve ID column if present; else create one
  orig_names <- names(df)
  names(df) <- clean_names_unique(orig_names)
  
  possible_id <- intersect(tolower(names(df)), c("row","sample","sample_id","id","case","case_id"))
  if (length(possible_id) == 0) {
    df$._row_id <- seq_len(nrow(df))
    id_col <- "._row_id"
  } else {
    id_col <- names(df)[match(possible_id[1], tolower(names(df)))]
  }
  
  # Column map
  name_map <- tibble::tibble(original_name = orig_names, clean_name = names(df))
  readr::write_csv(name_map, file.path(outdir, "00_column_name_map.csv"))
  
  # Drop identical duplicate columns created by make.unique
  dup_groups <- split(names(df), sub("__dup__.*$", "", names(df)))
  dup_groups <- dup_groups[vapply(dup_groups, length, integer(1)) > 1]
  
  dropped_cols <- character(0)
  qc_notes <- c()
  
  if (length(dup_groups) > 0) {
    for (g in names(dup_groups)) {
      cols <- dup_groups[[g]]
      base_col <- cols[1]
      identical_all <- all(vapply(cols[-1], function(cc) col_identical(df[[base_col]], df[[cc]]), logical(1)))
      if (identical_all) {
        drop_these <- cols[-1]
        df[drop_these] <- NULL
        dropped_cols <- c(dropped_cols, drop_these)
        qc_notes <- c(qc_notes, sprintf("Dropped identical duplicate columns for '%s': %s", g, paste(drop_these, collapse = ", ")))
      } else {
        qc_notes <- c(qc_notes, sprintf("Duplicate-name columns for '%s' exist but are NOT identical; kept all.", g))
      }
    }
  }
  
  n <- nrow(df); p <- ncol(df)
  
  # Missingness
  missing_by_var <- tibble::tibble(
    variable = names(df),
    n_missing = vapply(df, function(x) {
      if (is.character(x)) sum(is.na(x) | trimws(x) == "")
      else sum(is.na(x))
    }, integer(1)),
    pct_missing = 100 * n_missing / n
  ) %>% arrange(desc(pct_missing))
  
  readr::write_csv(missing_by_var, file.path(outdir, "02_missingness_by_variable.csv"))
  
  missing_by_row <- tibble::tibble(
    row_index = seq_len(n),
    n_missing = apply(df, 1, function(r) {
      rr <- as.list(r)
      sum(vapply(rr, function(x) {
        if (is.character(x)) is.na(x) || trimws(x) == ""
        else is.na(x)
      }, logical(1)))
    }),
    pct_missing = 100 * n_missing / p
  ) %>% arrange(desc(pct_missing))
  
  readr::write_csv(missing_by_row, file.path(outdir, "03_missingness_by_row.csv"))
  
  # Variable dictionary
  var_dict <- tibble::tibble(
    variable = names(df),
    role = vapply(names(df), infer_role, character(1)),
    inferred_type = vapply(df, infer_type, character(1)),
    n_unique_nonmissing = vapply(df, function(x) length(unique(x[!is.na(x)])), integer(1)),
    example_values = vapply(df, function(x) {
      vals <- unique(x[!is.na(x)])
      vals <- vals[1:min(5, length(vals))]
      paste(vals, collapse = "; ")
    }, character(1))
  ) %>%
    left_join(missing_by_var, by = "variable") %>%
    arrange(match(role, c("clinical_core","clinical_outcome","clinical_covariate","molecular_feature")),
            desc(pct_missing), variable)
  
  readr::write_csv(var_dict, file.path(outdir, "01_variable_dictionary.csv"))
  
  # Numeric summary
  numeric_vars <- var_dict %>%
    filter(inferred_type %in% c("numeric","integer_like","binary")) %>%
    pull(variable)
  
  summ_numeric <- function(x) {
    if (is.character(x)) x[trimws(x) == ""] <- NA
    xn <- suppressWarnings(as.numeric(x))
    xn <- xn[!is.na(xn)]
    if (length(xn) == 0) {
      return(tibble::tibble(
        n_nonmissing = 0, mean = NA_real_, sd = NA_real_, median = NA_real_,
        min = NA_real_, p01 = NA_real_, p25 = NA_real_, p75 = NA_real_,
        p99 = NA_real_, max = NA_real_, n_zeros = NA_integer_, n_negative = NA_integer_
      ))
    }
    qs <- stats::quantile(xn, probs = c(0.01,0.25,0.75,0.99), na.rm = TRUE, names = FALSE, type = 7)
    tibble::tibble(
      n_nonmissing = length(xn),
      mean   = mean(xn),
      sd     = stats::sd(xn),
      median = stats::median(xn),
      min    = min(xn),
      p01    = qs[1],
      p25    = qs[2],
      p75    = qs[3],
      p99    = qs[4],
      max    = max(xn),
      n_zeros    = sum(xn == 0),
      n_negative = sum(xn < 0)
    )
  }
  
  numeric_summary <- dplyr::bind_rows(lapply(numeric_vars, function(v) {
    out <- summ_numeric(df[[v]])
    out$variable <- v
    out
  })) %>%
    relocate(variable) %>%
    left_join(var_dict %>% select(variable, role, inferred_type, pct_missing), by = "variable") %>%
    arrange(match(role, c("clinical_core","clinical_outcome","clinical_covariate","molecular_feature")),
            desc(pct_missing), variable)
  
  readr::write_csv(numeric_summary, file.path(outdir, "04_numeric_summary.csv"))
  
  # Categorical summary
  cat_vars <- var_dict %>%
    filter(inferred_type %in% c("categorical","text","date_like")) %>%
    pull(variable)
  
  top_levels <- function(x, k = 10) {
    if (is.character(x)) x[trimws(x) == ""] <- NA
    x <- x[!is.na(x)]
    if (length(x) == 0) return("")
    tab <- sort(table(x), decreasing = TRUE)
    tab <- head(tab, k)
    paste(sprintf("%s (%d)", names(tab), as.integer(tab)), collapse = " | ")
  }
  
  categorical_summary <- tibble::tibble(
    variable = cat_vars,
    role = var_dict$role[match(cat_vars, var_dict$variable)],
    inferred_type = var_dict$inferred_type[match(cat_vars, var_dict$variable)],
    n_unique_nonmissing = var_dict$n_unique_nonmissing[match(cat_vars, var_dict$variable)],
    pct_missing = var_dict$pct_missing[match(cat_vars, var_dict$variable)],
    top_levels = vapply(cat_vars, function(v) top_levels(df[[v]], k = 10), character(1))
  ) %>%
    arrange(match(role, c("clinical_core","clinical_outcome","clinical_covariate","molecular_feature")),
            desc(pct_missing), variable)
  
  readr::write_csv(categorical_summary, file.path(outdir, "05_categorical_summary.csv"))
  
  # QC flags
  flag_lines <- c()
  add_flag <- function(...) flag_lines <<- c(flag_lines, sprintf(...))
  
  core_expected <- c("age_at_dx","age","race","stage","os_years","os_status","usc73")
  core_present <- intersect(core_expected, names(df))
  core_missing <- setdiff(core_expected, names(df))
  add_flag("Rows: %d | Columns: %d", nrow(df), ncol(df))
  if (length(core_present) > 0) add_flag("Clinical core present: %s", paste(core_present, collapse = ", "))
  if (length(core_missing) > 0) add_flag("Clinical core NOT found (check naming): %s", paste(core_missing, collapse = ", "))
  
  if ("os_years" %in% names(df)) {
    os <- suppressWarnings(as.numeric(df$os_years))
    add_flag("OS years: non-missing=%d; min=%.3f; max=%.3f; negatives=%d",
             sum(!is.na(os)), suppressWarnings(min(os, na.rm = TRUE)),
             suppressWarnings(max(os, na.rm = TRUE)), sum(os < 0, na.rm = TRUE))
  }
  if ("os_status" %in% names(df)) {
    st <- as.character(df$os_status)
    st[trimws(st) == ""] <- NA
    uniq <- unique(st[!is.na(st)])
    add_flag("OS status unique values (up to 20): %s", paste(head(uniq, 20), collapse = ", "))
  }
  if ("stage" %in% names(df)) {
    s <- as.character(df$stage); s[trimws(s) == ""] <- NA
    uniqs <- unique(s[!is.na(s)])
    add_flag("Stage unique values (up to 30): %s", paste(head(uniqs, 30), collapse = ", "))
  }
  
  const_cols <- names(df)[vapply(df, function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))]
  if (length(const_cols) > 0) add_flag("Constant columns (<=1 unique non-missing): %s", paste(const_cols, collapse = ", "))
  
  high_miss <- missing_by_var %>% filter(pct_missing >= 30) %>% pull(variable)
  if (length(high_miss) > 0) add_flag("High-missing columns (>=30%%): %d (see 02_missingness_by_variable.csv)", length(high_miss))
  
  if (id_col %in% names(df)) {
    if (any(duplicated(df[[id_col]]))) add_flag("WARNING: Duplicate IDs detected in '%s'.", id_col)
  }
  
  if (length(dropped_cols) > 0) add_flag("Dropped identical duplicate columns: %s", paste(dropped_cols, collapse = ", "))
  if (length(qc_notes) > 0) add_flag("Duplicate-column notes:\n- %s", paste(qc_notes, collapse = "\n- "))
  
  writeLines(c("=== DATA QC FLAGS / NOTES ===", flag_lines), con = file.path(outdir, "06_qc_flags.txt"))
  
  # Optional plots
  if (plots) {
    top_miss <- missing_by_var %>% slice_head(n = min(40, nrow(missing_by_var)))
    p1 <- ggplot(top_miss, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
      geom_col() + coord_flip() +
      labs(title = "Top missingness by variable", x = "Variable", y = "% missing")
    ggsave(file.path(outdir, "plots", "missingness_top40.png"), p1, width = 9, height = 7, dpi = 200)
    
    key_num <- intersect(c("usc73","age","age_at_dx","os_years"), names(df))
    for (v in key_num) {
      x <- suppressWarnings(as.numeric(df[[v]]))
      dd <- tibble::tibble(value = x)
      p <- ggplot(dd, aes(x = value)) +
        geom_histogram(bins = 30) +
        labs(title = paste("Histogram:", v), x = v, y = "Count")
      ggsave(file.path(outdir, "plots", paste0("hist_", v, ".png")), p, width = 7, height = 5, dpi = 200)
    }
  }
  
  invisible(list(df = df, var_dict = var_dict, missing_by_var = missing_by_var))
}

# -------------------------------
# Read sheets + auto-select
# -------------------------------
sheets <- readxl::excel_sheets(INPUT)
if (length(sheets) == 0) stop("No sheets found in Excel file.", call. = FALSE)

message("Found sheets:\n - ", paste(sheets, collapse = "\n - "))

read_sheet_safe <- function(sh) {
  tryCatch(
    readxl::read_excel(INPUT, sheet = sh, guess_max = GUESSMAX),
    error = function(e) NULL
  )
}

if (!is.na(SHEET_IN)) {
  if (!(SHEET_IN %in% sheets)) stop("Requested sheet not found: ", SHEET_IN, call. = FALSE)
  chosen_sheet <- SHEET_IN
  message("Using user-specified sheet: ", chosen_sheet)
} else {
  # Score each sheet
  sheet_info <- purrr::map_dfr(sheets, function(sh) {
    df <- read_sheet_safe(sh)
    if (is.null(df)) {
      tibble::tibble(sheet = sh, clin_hits = NA_integer_, ncol = NA_integer_, nrow = NA_integer_, score = NA_real_)
    } else {
      si <- score_sheet(df)
      tibble::tibble(sheet = sh) %>% bind_cols(si)
    }
  }) %>% arrange(desc(score))
  
  readr::write_csv(sheet_info, file.path(OUTDIR, "sheet_manifest_scores.csv"))
  
  chosen_sheet <- sheet_info$sheet[which.max(sheet_info$score)]
  message("\nAuto-selected sheet: ", chosen_sheet)
  message("Wrote sheet scoring to: ", file.path(OUTDIR, "sheet_manifest_scores.csv"))
}

df <- read_sheet_safe(chosen_sheet)
if (is.null(df)) stop("Failed to read chosen sheet: ", chosen_sheet, call. = FALSE)

# -------------------------------
# Profile chosen sheet
# -------------------------------
sheet_outdir <- file.path(OUTDIR, paste0("sheet__", janitor::make_clean_names(chosen_sheet)))
res <- profile_df(df, sheet_outdir, plots = PLOTS)

message("\nDONE.")
message("Outputs written under: ", normalizePath(sheet_outdir))
message("Start with:")
message(" - 06_qc_flags.txt")
message(" - 01_variable_dictionary.csv")
message(" - 02_missingness_by_variable.csv")
message(" - 04_numeric_summary.csv")
