#' =============================================================================
#' Compositional Analysis of a Historical Landscape-Painting Corpus
#' -----------------------------------------------------------------------------
#' Research question
#'   Does the depicted land-cover composition of European landscape paintings
#'   change across 1611-1899 (i.e. is there a shifting-baseline signal)?
#'
#' Design principle
#'   Each artwork is ONE compositional observation. The composition is the
#'   response. Only variables that are NOT functions of the composition may
#'   appear on the right-hand side -- hence no DAFOR score and no pixel-derived
#'   habitat class as predictors.
#'
#' Pipeline stages
#'   STAGE 1 (sections 2-10) runs in seconds and requires no Stan. Run it first:
#'           it reproduces every descriptive and frequentist number in the paper.
#'   STAGE 2 (section 11) fits the Bayesian multivariate models with brms and is
#'           the expensive part. Toggle with the RUN_STAGE2 flag / env var.
#'
#' Inputs (read from DATA_DIR; see section 1)
#'   composition_counts.csv          per-artwork integer pixel counts per class
#'   clean_metadata.xlsx             Image_Name, median_date, latitude,
#'                                   longitude, painter_name,
#'                                   painter_school_tradition and
#'                                   dominant_landscape_type (the EUNIS-aligned
#'                                   scene type assigned independently of the
#'                                   segmentation model)
#'   R_pooled_rownorm.csv            5x5 row-normalised confusion matrix,
#'                                   no header, order = wooded, open_veg,
#'                                   substrate, water, human
#'   per_painting_confusion_long.csv stem, true_class, pred_class, n
#'
#' Outputs (written to OUT_DIR; see sections 14-15)
#'   analysis_frame.csv, diagnostic_battery.csv, diagnostic_battery_nnls.csv,
#'   within_scene_type.csv, variance_partition.csv, bootstrap_quantiles.csv,
#'   school_mapping.csv, scene_type_counts.csv, works_per_painter.csv,
#'   posterior_summaries.csv, loo_comparison.txt, session_info.txt,
#'   run_manifest.csv, ppc_*.png
#'
#' Reproducibility
#'   Deterministic given the same inputs and RNG_SEED.
#'
#' Author  : Ahmet Selim Esmer
#' Licence : 
#' Style   : tidyverse style guide (https://style.tidyverse.org)
#' =============================================================================


# =============================================================================
# SECTION 1 - CONFIGURATION, DEPENDENCIES AND HELPERS
# -----------------------------------------------------------------------------
# Purpose : Declare every dependency, resolve all filesystem locations from the
#           environment, fix the RNG seed, and define the logging and assertion
#           helpers used by the inter-section test blocks.
# Inputs  : Environment variables PAINTINGS_DATA, PAINTINGS_OUT, PAINTINGS_STAGE2.
# Outputs : Global configuration constants, OUT_DIR created on disk, helper
#           functions log_step(), check_that() and expect_cols().
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(stringi)
  library(tidyr)
  library(tibble)
  library(nnls)
  library(brms)
  library(lme4)
})

## --- Runtime flags -----------------------------------------------------------
# Stage 2 (brms) is opt-out rather than opt-in so a bare `Rscript` run is cheap.
#' Parse a truthy environment variable ("TRUE"/"true"/"1"/"yes" all count).
parse_flag <- function(name, default = "FALSE") {
  value <- tolower(trimws(Sys.getenv(name, unset = default)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n", "")) return(FALSE)
  stop(sprintf("environment variable %s must be boolean, got '%s'", name, value),
       call. = FALSE)
}

RUN_STAGE2           <- parse_flag("PAINTINGS_STAGE2")
RUN_PIXEL_COMPARISON <- TRUE   # evaluate the legacy pixel rule against the expert scheme
RUN_PRIOR_PREDICTIVE <- TRUE   # refits fit_A with sample_prior = "only" (slow)

## --- Filesystem --------------------------------------------------------------
# Repo-relative defaults; override with environment variables in CI or locally.
DATA_DIR <- Sys.getenv("PAINTINGS_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("PAINTINGS_OUT",  unset = "outputs")
FIG_DIR  <- file.path(OUT_DIR, "figures")


DATA_DIR <- Sys.getenv("PAINTINGS_DATA", unset = "C:/Users/186Poplar/Documents/thesis_data")
OUT_DIR  <- Sys.getenv("PAINTINGS_OUT",  unset = "C:/Users/186Poplar/Documents/thesis_data/outputs")

FILE_COUNTS    <- "composition_counts.csv"
FILE_METADATA  <- "clean_metadata.xlsx"
FILE_CONFUSION <- "R_pooled_rownorm.csv"
FILE_PER_PAINT <- "per_painting_confusion_long.csv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## --- Reproducibility ---------------------------------------------------------
RNG_SEED <- 42L
set.seed(RNG_SEED)

## --- Analysis constants ------------------------------------------------------
# Expected corpus shape. These are assertions, not documentation: if the inputs
# change, the pipeline must fail loudly here rather than silently downstream.
EXPECTED_N              <- 238L
EXPECTED_SCHOOL_LEVELS  <- 10L
EXPECTED_DUTCH_NORTHERN <- 115L
EXPECTED_DUTCH_CORE     <- 110L
DATE_MIN                <- 1600
DATE_MAX                <- 1900
CENTURY_BREAKS          <- c(1600, 1700, 1800, 1900)
CENTURY_LABELS          <- c("17th", "18th", "19th")

DETECTION_LIMIT_PX <- 0.65   # sub-pixel detection limit for zero replacement
RIDGE_LAMBDA       <- 50     # closure weight in the constrained inversion
RIDGE_RHO          <- 0.30   # shrinkage toward the observed composition
N_BOOTSTRAP        <- 300L   # painting-level bootstrap draws of R
POWER_NSIM         <- 400L   # simulations per effect size in the design analysis
POWER_EFFECTS      <- seq(0, 0.5, by = 0.05)
MIN_SUBSET_N       <- 20L    # smallest subset the zero-exclusion check will fit
RESPONSES          <- c("b1", "b2", "b3", "b4")

## --- Package version guard ---------------------------------------------------
REQUIRED_PACKAGES <- c("dplyr", "readr", "readxl", "stringr", "stringi",
                       "tidyr", "tibble", "nnls", "lme4")
STAGE2_PACKAGES   <- c("brms", "loo", "posterior", "ggplot2")

assert_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: install.packages(c(",
         paste0('"', missing_pkgs, '"', collapse = ", "), "))",
         call. = FALSE)
  }
  invisible(TRUE)
}
assert_packages(REQUIRED_PACKAGES)
if (RUN_STAGE2) assert_packages(STAGE2_PACKAGES)

## --- Logging and assertion helpers -------------------------------------------
RUN_STARTED_AT <- Sys.time()

log_step <- function(...) {
  elapsed <- difftime(Sys.time(), RUN_STARTED_AT, units = "secs")
  message(sprintf("[%6.1fs] %s", as.numeric(elapsed), paste0(...)))
  invisible(NULL)
}

#' Assert a condition and fail with a section-scoped message.
check_that <- function(condition, message_text) {
  if (!isTRUE(all(condition, na.rm = FALSE))) {
    stop("VALIDATION FAILED: ", message_text, call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a data frame contains the named columns.
expect_cols <- function(data, cols, object_name = deparse(substitute(data))) {
  missing_cols <- setdiff(cols, names(data))
  check_that(length(missing_cols) == 0L,
             sprintf("%s is missing column(s): %s",
                     object_name, paste(missing_cols, collapse = ", ")))
}

#' Assert that a numeric object is entirely finite (no NA, NaN or Inf).
expect_finite <- function(x, object_name = deparse(substitute(x))) {
  check_that(all(is.finite(as.matrix(x))),
             sprintf("%s contains non-finite values (NA, NaN or Inf)", object_name))
}

resolve_input <- function(filename) {
  path <- file.path(DATA_DIR, filename)
  check_that(file.exists(path),
             sprintf("input file not found: %s (set PAINTINGS_DATA)", path))
  path
}

as_century <- function(year) {
  cut(year, CENTURY_BREAKS, labels = CENTURY_LABELS, right = FALSE)
}

log_step("configuration loaded | data = '", DATA_DIR, "' | out = '", OUT_DIR,
         "' | stage2 = ", RUN_STAGE2)


# ---- TEST BLOCK 1 -> validates SECTION 1 ------------------------------------
local({
  check_that(dir.exists(DATA_DIR), "DATA_DIR does not exist")
  check_that(dir.exists(OUT_DIR),  "OUT_DIR was not created")
  check_that(dir.exists(FIG_DIR),  "FIG_DIR was not created")
  check_that(is.logical(RUN_STAGE2) && !is.na(RUN_STAGE2),
             "RUN_STAGE2 must resolve to TRUE or FALSE")
  check_that(EXPECTED_N > 0L && EXPECTED_SCHOOL_LEVELS > 0L,
             "corpus constants must be positive")
  check_that(RIDGE_LAMBDA > 0 && RIDGE_RHO > 0,
             "regularisation constants must be positive")
  check_that(is.function(check_that) && is.function(log_step),
             "helper functions were not defined")
  log_step("TEST 1 passed - configuration and helpers valid")
})


# =============================================================================
# SECTION 2 - DATA INGEST, HARMONISATION AND THE ANALYSIS FRAME
# -----------------------------------------------------------------------------
# Purpose : Join pixel counts to curatorial metadata, apply the inclusion rule,
#           derive standardised covariates and the canonical school assignment.
# Inputs  : composition_counts.csv, clean_metadata.xlsx.
# Outputs : `analysis_data` (EXPECTED_N rows) and analysis_frame.csv, which is
#           the single source of truth consumed by all downstream figure code.
# =============================================================================

#' Canonical painter-school assignment (the random-effect grouping factor).
#'
#' Diacritics are stripped before matching. Matching the literal "dusseldorf"
#' against a UTF-8 "d??sseldorf" silently fails under some locales, which
#' previously routed three artworks into "Unspecified / Other" and desynchronised
#' the published figure from the fitted models.
school_clean <- function(x) {
  s <- str_trim(str_to_lower(str_split_fixed(x, ",", 2)[, 1]))
  s <- stri_trans_general(s, "Latin-ASCII")
  case_when(
    str_detect(s, "dutch|hague|dusseldorf|danish")      ~ "Dutch & Northern",
    str_detect(s, "realism|barbizon")                   ~ "Realism & Barbizon",
    str_detect(s, "classicism|classical|neoclassicism") ~ "Classicism & Neoclassicism",
    str_detect(s, "romanticism")                        ~ "Romanticism",
    str_detect(s, "impressionism")                      ~ "Impressionism",
    str_detect(s, "rococo|arcad")                       ~ "Rococo & Arcadian",
    str_detect(s, "italianate|vedut|venetian")          ~ "Italianate & Vedute",
    str_detect(s, "baroque")                            ~ "Baroque",
    str_detect(s, "spanish")                            ~ "Spanish School",
    TRUE                                                ~ "Unspecified / Other"
  )
}

#' Narrow within-tradition membership, used only for the single-tradition refit.
#'
#' Deliberately stricter than `school_clean()`: Calame was Swiss and the
#' Duesseldorf academy German, so those works share a northern convention (and
#' therefore a random-effect level) without belonging to the Dutch Golden Age
#' tradition. This is why the restricted model has EXPECTED_DUTCH_CORE rows.
is_dutch_core <- function(x) {
  s <- stri_trans_general(str_to_lower(str_split_fixed(x, ",", 2)[, 1]),
                          "Latin-ASCII")
  str_detect(s, "dutch|hague")
}

counts_raw <- read_csv(resolve_input(FILE_COUNTS), show_col_types = FALSE)
metadata_raw <- read_excel(resolve_input(FILE_METADATA)) |>
  rename(habitat_expert = dominant_landscape_type)

expect_cols(counts_raw, "Image_Name")
expect_cols(metadata_raw,
            c("Image_Name", "median_date", "latitude", "longitude",
              "painter_name", "painter_school_tradition", "habitat_expert"))
check_that(anyDuplicated(counts_raw$Image_Name) == 0L,
           "duplicate Image_Name in composition_counts.csv")
check_that(anyDuplicated(metadata_raw$Image_Name) == 0L,
           "duplicate Image_Name in the metadata workbook")

joined_raw <- inner_join(counts_raw, metadata_raw, by = "Image_Name")
log_step(sprintf("joined %d of %d artworks", nrow(joined_raw), nrow(counts_raw)))

analysis_data <- joined_raw |>
  mutate(
    median_date = as.numeric(median_date),
    latitude    = as.numeric(latitude),
    longitude   = as.numeric(longitude)
  ) |>
  filter(
    !is.na(median_date), median_date >= DATE_MIN, median_date < DATE_MAX,
    !is.na(latitude), !is.na(longitude)
  ) |>
  mutate(
    date_sc    = as.numeric(scale(median_date)),
    lat_sc     = as.numeric(scale(latitude)),
    lon_sc     = as.numeric(scale(longitude)),
    painter    = factor(painter_name),
    school     = factor(school_clean(painter_school_tradition)),
    dutch_core = is_dutch_core(painter_school_tradition),
    habitat    = factor(habitat_expert)
  )

# The frame every figure script reads. Written here, immediately after the
# analysis set is defined, so no downstream consumer can recompute and diverge.
analysis_frame <- analysis_data |>
  select(Image_Name, median_date, latitude, longitude,
         painter, school, dutch_core, habitat)
write_csv(analysis_frame, file.path(OUT_DIR, "analysis_frame.csv"))

# Auditable raw-label -> assigned-category mapping (reported in the appendix).
school_mapping <- analysis_data |>
  mutate(raw_label = stri_trans_general(
    str_to_lower(str_split_fixed(painter_school_tradition, ",", 2)[, 1]),
    "Latin-ASCII"
  )) |>
  count(school, raw_label, name = "n") |>
  arrange(school, desc(n))
write_csv(school_mapping, file.path(OUT_DIR, "school_mapping.csv"))

log_step(sprintf("analysis set: n = %d, %d-%d, %d painters, %d schools",
                 nrow(analysis_data), min(analysis_data$median_date),
                 max(analysis_data$median_date),
                 nlevels(analysis_data$painter), nlevels(analysis_data$school)))


# ---- TEST BLOCK 2 -> validates SECTION 2 ------------------------------------
local({
  check_that(nrow(analysis_data) == EXPECTED_N,
             sprintf("analysis set has %d rows, expected %d",
                     nrow(analysis_data), EXPECTED_N))
  check_that(nlevels(analysis_data$school) == EXPECTED_SCHOOL_LEVELS,
             sprintf("school has %d levels, expected %d",
                     nlevels(analysis_data$school), EXPECTED_SCHOOL_LEVELS))
  check_that(sum(analysis_data$school == "Dutch & Northern") == EXPECTED_DUTCH_NORTHERN,
             "Dutch & Northern group size drifted (check the diacritic handling)")
  check_that(sum(analysis_data$dutch_core) == EXPECTED_DUTCH_CORE,
             "dutch_core subset size drifted")
  check_that(anyDuplicated(analysis_data$Image_Name) == 0L,
             "duplicate artworks survived the join")
  
  expect_finite(analysis_data[c("median_date", "latitude", "longitude",
                                "date_sc", "lat_sc", "lon_sc")])
  check_that(!anyNA(analysis_data$school) && !anyNA(analysis_data$habitat),
             "NA present in school or habitat")
  check_that(abs(mean(analysis_data$date_sc)) < 1e-8 &&
               abs(sd(analysis_data$date_sc) - 1) < 1e-8,
             "date_sc is not standardised to mean 0, sd 1")
  
  century_n <- as.integer(table(as_century(analysis_data$median_date)))
  check_that(sum(century_n) == EXPECTED_N, "century bins do not partition the corpus")
  
  check_that(file.exists(file.path(OUT_DIR, "analysis_frame.csv")),
             "analysis_frame.csv was not written")
  
  # A scene-type level with very few artworks yields an unidentified fixed
  # effect. Warn rather than fail: the correct action is editorial, not runtime.
  sparse_levels <- names(which(table(analysis_data$habitat) < 5L))
  if (length(sparse_levels) > 0L) {
    warning("scene types with n < 5 (coefficients will be weakly identified): ",
            paste(sparse_levels, collapse = ", "), call. = FALSE)
  }
  log_step("TEST 2 passed - analysis frame valid (n = ", nrow(analysis_data), ")")
})


# =============================================================================
# SECTION 3 - POOLING RAW CLASSES INTO FIVE ECOLOGICAL PARTS
# -----------------------------------------------------------------------------
# Purpose : Aggregate the 15 segmentation classes into five parts chosen on both
#           ecological meaning and measurement reliability.
# Inputs  : `analysis_data` pixel-count columns (prefix "n_").
# Outputs : `part_counts` (n x 5 integer matrix), `part_totals`, `P_observed`.
# =============================================================================

# tree_conical (2.6% IoU, absent from 96.6% of artworks) and path_road (2.2%
# IoU) are absorbed into the aggregates they belong to rather than kept as
# separate parts: 74% of true tree_conical pixels are predicted tree_broadleaf,
# a confusion that is within-pool and therefore correct at the ecological level.
PARTS <- list(
  wooded    = c("wooded_mass", "tree_broadleaf", "tree_conical"),
  open_veg  = c("grass", "shrub_bush"),
  substrate = c("earth", "rock", "mountain"),
  water     = c("Water"),
  human     = c("building", "path_road", "person")
)
PART_NAMES <- names(PARTS)

expect_cols(analysis_data, paste0("n_", unlist(PARTS, use.names = FALSE)))

part_counts <- vapply(
  PARTS,
  function(classes) {
    rowSums(as.matrix(as.data.frame(analysis_data)[, paste0("n_", classes),
                                                   drop = FALSE]))
  },
  numeric(nrow(analysis_data))
)
colnames(part_counts) <- PART_NAMES
part_totals <- rowSums(part_counts)
P_observed  <- part_counts / part_totals

log_step("mean cover: ",
         paste(sprintf("%s %.1f%%", PART_NAMES, 100 * colMeans(P_observed)),
               collapse = "  "))


# ---- TEST BLOCK 3 -> validates SECTION 3 ------------------------------------
local({
  check_that(identical(dim(part_counts), c(EXPECTED_N, length(PARTS))),
             "part_counts has unexpected dimensions")
  check_that(identical(colnames(P_observed), PART_NAMES),
             "part column names are out of order")
  expect_finite(part_counts)
  expect_finite(P_observed)
  check_that(all(part_counts >= 0), "negative pixel counts present")
  check_that(all(part_totals > 0), "artwork(s) with zero terrestrial pixels")
  check_that(max(abs(rowSums(P_observed) - 1)) < 1e-9,
             "compositions do not close to 1")
  check_that(all(P_observed >= 0 & P_observed <= 1),
             "compositional shares outside [0, 1]")
  log_step("TEST 3 passed - five-part composition closes to 1")
})


# =============================================================================
# SECTION 4 - LEGACY PIXEL-RULE HABITAT (EVALUATED, NEVER USED AS A PREDICTOR)
# -----------------------------------------------------------------------------
# Purpose : Recompute the threshold habitat rule solely so it can be scored
#           against the independent expert classification.
# Inputs  : `P_observed`, `analysis_data$habitat`.
# Outputs : `analysis_data$hab_pixel`, `pixel_rule_agreement` (list).
# =============================================================================

# The rule is a deterministic function of the response, which is precisely what
# made the original model circular. It is retained here as a classifier under
# evaluation and is never admitted to the right-hand side of any model.
EXPERT_TO_PIXEL_SCHEME <- c(
  "woodland landscape"            = "Woodland & Forest",
  "open countryside"              = "Open Semi-Natural",
  "river/lakeside landscape"      = "Aquatic & Riparian",
  "coastal landscape"             = "Aquatic & Riparian",
  "mountain/rocky landscape"      = "Sparsely Vegetated / Montane",
  "urban/architectural landscape" = "Anthropogenic / Cultural"
)

classify_pixel_rule <- function(P) {
  shares <- as.data.frame(P)
  factor(with(shares, case_when(
    water     > 0.25     ~ "Aquatic & Riparian",
    wooded    > 0.50     ~ "Woodland & Forest",
    human     > 0.15     ~ "Anthropogenic / Cultural",
    substrate > open_veg ~ "Sparsely Vegetated / Montane",
    TRUE                 ~ "Open Semi-Natural"
  )))
}

#' Cohen's kappa and raw agreement for two aligned categorical vectors.
cohens_kappa <- function(x, y) {
  usable <- !is.na(x) & !is.na(y)
  observed <- table(x[usable], y[usable])
  p_observed <- sum(diag(observed)) / sum(observed)
  p_expected <- sum(rowSums(observed) * colSums(observed)) / sum(observed)^2
  list(agreement = p_observed,
       kappa     = (p_observed - p_expected) / (1 - p_expected),
       n         = sum(observed),
       table     = observed)
}

pixel_rule_agreement <- NULL
if (RUN_PIXEL_COMPARISON) {
  analysis_data$hab_pixel <- classify_pixel_rule(P_observed)
  pixel_rule_agreement <- cohens_kappa(
    EXPERT_TO_PIXEL_SCHEME[as.character(analysis_data$habitat)],
    as.character(analysis_data$hab_pixel)
  )
  log_step(sprintf("expert vs pixel rule: agreement %.1f%%, kappa %.3f (n = %d)",
                   100 * pixel_rule_agreement$agreement,
                   pixel_rule_agreement$kappa, pixel_rule_agreement$n))
}

scene_type_counts <- analysis_data |>
  count(habitat, name = "n") |>
  arrange(desc(n))
write_csv(scene_type_counts, file.path(OUT_DIR, "scene_type_counts.csv"))


# ---- TEST BLOCK 4 -> validates SECTION 4 ------------------------------------
local({
  check_that(sum(scene_type_counts$n) == EXPECTED_N,
             "scene-type counts do not sum to the corpus size")
  if (RUN_PIXEL_COMPARISON) {
    check_that("hab_pixel" %in% names(analysis_data),
               "hab_pixel was not attached")
    check_that(!anyNA(analysis_data$hab_pixel),
               "pixel rule produced NA classifications")
    check_that(is.finite(pixel_rule_agreement$kappa),
               "kappa for the pixel rule is not finite")
    check_that(pixel_rule_agreement$kappa >= -1 && pixel_rule_agreement$kappa <= 1,
               "kappa outside the valid [-1, 1] range")
    check_that(pixel_rule_agreement$n > 0, "no artworks matched between schemes")
  }
  log_step("TEST 4 passed - pixel-rule evaluation complete")
})


# =============================================================================
# SECTION 5 - CONFUSION-MATRIX BIAS CORRECTION
# -----------------------------------------------------------------------------
# Purpose : Invert the row-normalised confusion matrix under simplex constraints
#           to obtain bias-corrected compositions (Card 1982; Olofsson 2014).
# Inputs  : `P_observed`, R_pooled_rownorm.csv.
# Outputs : `R_matrix`, `P_ridge` (primary correction), `P_nnls` (failure demo).
# =============================================================================

#' Constrained inversion of the confusion matrix.
#'
#' Plain NNLS returns SPARSE solutions: it places mass on the simplex boundary
#' by construction, converting small non-zero parts into exact zeros. Those
#' zeros are then zero-replaced and land at the extremes of the log-ratio scale,
#' which manufactured a spurious woodland decline that survived both bootstrap
#' resampling and detection-limit sensitivity. The ridge term shrinks toward the
#' observed composition and keeps solutions off the boundary. Always report
#' "none" alongside "ridge"; "nnls" exists to demonstrate the failure mode.
correct_composition <- function(P, R, method = c("ridge", "nnls", "none"),
                                lambda = RIDGE_LAMBDA, rho = RIDGE_RHO) {
  method <- match.arg(method)
  if (method == "none") return(P)
  
  k <- ncol(R)
  design <- rbind(t(R), sqrt(lambda) * rep(1, k))
  if (method == "ridge") design <- rbind(design, sqrt(rho) * diag(k))
  
  corrected <- t(apply(P, 1, function(p) {
    target <- c(p, sqrt(lambda))
    if (method == "ridge") target <- c(target, sqrt(rho) * p)
    solution <- nnls(design, target)$x
    if (sum(solution) > 0) solution / sum(solution) else p
  }))
  colnames(corrected) <- colnames(P)
  corrected
}

R_matrix <- as.matrix(read.csv(resolve_input(FILE_CONFUSION), header = FALSE))
dimnames(R_matrix) <- list(PART_NAMES, PART_NAMES)

P_ridge <- correct_composition(P_observed, R_matrix, "ridge")
P_nnls  <- correct_composition(P_observed, R_matrix, "nnls")

log_step(sprintf("zero cells -- observed %d | ridge %d | nnls %d",
                 sum(P_observed <= 0), sum(P_ridge <= 1e-9),
                 sum(P_nnls <= 1e-9)))


# ---- TEST BLOCK 5 -> validates SECTION 5 ------------------------------------
local({
  check_that(identical(dim(R_matrix), c(5L, 5L)),
             "confusion matrix is not 5x5")
  expect_finite(R_matrix)
  check_that(all(R_matrix >= 0), "negative entries in the confusion matrix")
  check_that(max(abs(rowSums(R_matrix) - 1)) < 1e-6,
             "confusion matrix rows are not normalised to 1")
  # base::kappa() is the matrix condition number here, not Cohen's kappa.
  check_that(kappa(t(R_matrix), exact = TRUE) < 30,
             "confusion matrix is ill-conditioned; inversion is unstable")
  
  for (nm in c("P_ridge", "P_nnls")) {
    mat <- get(nm)
    check_that(identical(dim(mat), dim(P_observed)),
               sprintf("%s has unexpected dimensions", nm))
    expect_finite(mat)
    check_that(all(mat >= 0), sprintf("%s contains negative shares", nm))
    check_that(max(abs(rowSums(mat) - 1)) < 1e-9,
               sprintf("%s rows do not close to 1", nm))
  }
  check_that(sum(P_nnls <= 1e-9) >= sum(P_ridge <= 1e-9),
             "ridge produced more boundary zeros than NNLS; check regularisation")
  log_step("TEST 5 passed - corrections stay on the simplex")
})


# =============================================================================
# SECTION 6 - ZERO REPLACEMENT AND ISOMETRIC LOG-RATIO BALANCES
# -----------------------------------------------------------------------------
# Purpose : Replace structural zeros multiplicatively, then map each composition
#           to four orthonormal ilr balances from an a priori binary partition.
# Inputs  : A compositional matrix and its per-artwork pixel totals.
# Outputs : `make_balances()` returning a data frame of b1-b4.
# =============================================================================

# Sequential binary partition fixed a priori, each balance answering one question:
#   b1 natural vs human      (wooded, open_veg, substrate, water) | (human)
#   b2 vegetation vs abiotic (wooded, open_veg)                   | (substrate, water)
#   b3 wooded vs open        (wooded)                             | (open_veg)
#   b4 substrate vs water    (substrate)                          | (water)
# Orthonormal and subcompositionally coherent (Aitchison 1986; Egozcue 2003).
BALANCE_PARTITION <- list(
  b1 = list(num = c("wooded", "open_veg", "substrate", "water"), den = "human"),
  b2 = list(num = c("wooded", "open_veg"), den = c("substrate", "water")),
  b3 = list(num = "wooded",    den = "open_veg"),
  b4 = list(num = "substrate", den = "water")
)

# Parts entering each balance; used by the zero-exclusion diagnostic.
BALANCE_PARTS <- lapply(BALANCE_PARTITION, function(p) c(p$num, p$den))

#' Multiplicative replacement of structural zeros at the per-artwork detection limit.
replace_zeros <- function(P, totals, px = DETECTION_LIMIT_PX) {
  detection_limit <- px / totals
  is_zero <- P <= 0
  limits <- matrix(rep(detection_limit, ncol(P)), nrow = nrow(P))
  replaced <- P * (1 - rowSums(is_zero * limits))
  replaced[is_zero] <- limits[is_zero]
  replaced / rowSums(replaced)
}

#' One orthonormal ilr balance between two groups of parts.
balance <- function(P, numerator, denominator) {
  r <- length(numerator)
  s <- length(denominator)
  geometric_mean <- function(M) exp(rowMeans(log(M)))
  sqrt(r * s / (r + s)) *
    log(geometric_mean(P[, numerator, drop = FALSE]) /
          geometric_mean(P[, denominator, drop = FALSE]))
}

#' All four balances for a composition matrix, after zero replacement.
make_balances <- function(P, totals, px = DETECTION_LIMIT_PX) {
  replaced <- replace_zeros(P, totals, px)
  as.data.frame(lapply(BALANCE_PARTITION, function(p) {
    balance(replaced, p$num, p$den)
  }))
}

balances_observed <- make_balances(P_observed, part_totals)


# ---- TEST BLOCK 6 -> validates SECTION 6 ------------------------------------
local({
  replaced <- replace_zeros(P_observed, part_totals)
  check_that(all(replaced > 0), "zero replacement left non-positive shares")
  check_that(max(abs(rowSums(replaced) - 1)) < 1e-9,
             "zero-replaced composition does not close to 1")
  
  check_that(identical(names(balances_observed), RESPONSES),
             "balance names do not match RESPONSES")
  check_that(nrow(balances_observed) == EXPECTED_N,
             "balances have the wrong number of rows")
  expect_finite(balances_observed)
  
  # Subcompositional coherence: b3 must be identical whether computed from the
  # full five-part composition or from the wooded/open subcomposition alone.
  # Checked on artworks with no structural zeros, where the identity is exact.
  zero_free <- rowSums(P_observed[, c("wooded", "open_veg"), drop = FALSE] <= 0) == 0
  check_that(sum(zero_free) > 0L, "no zero-free artworks available for the check")
  subcomp <- P_observed[zero_free, c("wooded", "open_veg"), drop = FALSE]
  subcomp <- subcomp / rowSums(subcomp)
  b3_sub <- balance(subcomp, "wooded", "open_veg")
  check_that(max(abs(b3_sub - balances_observed$b3[zero_free])) < 1e-8,
             "b3 is not subcompositionally coherent")
  log_step("TEST 6 passed - balances finite and coherent")
})


# =============================================================================
# SECTION 7 - CLUSTER-ROBUST OLS UTILITY
# -----------------------------------------------------------------------------
# Purpose : Estimate one coefficient with painter-clustered (sandwich) standard
#           errors, degrading gracefully on rank-deficient subsets.
# Inputs  : A design matrix, a response vector and a clustering vector.
# Outputs : `ols_cluster()` returning est, se, rank and a status note.
# =============================================================================

# The painter, not the painting, is the unit of replication: a single artist
# contributing 16 works does not supply 16 independent observations.
ols_cluster <- function(design, y, cluster, coef_name = "date") {
  degenerate <- list(est = NA_real_, se = NA_real_, rank = 0L,
                     note = "degenerate")
  
  design <- as.matrix(design)
  if (is.null(colnames(design))) {
    colnames(design) <- paste0("v", seq_len(ncol(design)))
  }
  
  # Drop rows carrying any non-finite value.
  usable <- is.finite(y) & apply(design, 1, function(row) all(is.finite(row)))
  if (sum(usable) < 5L) return(degenerate)
  design  <- design[usable, , drop = FALSE]
  y       <- y[usable]
  cluster <- cluster[usable]
  
  # Drop constant columns (NA-safe); the intercept is protected by name.
  column_var <- apply(design, 2, stats::var)
  keep <- !is.na(column_var) & column_var > 0
  keep[colnames(design) == "(int)"] <- TRUE
  if (!any(keep)) return(degenerate)
  design <- design[, keep, drop = FALSE]
  
  # Drop linearly dependent columns via QR pivoting.
  decomposition <- qr(design)
  if (decomposition$rank < 1L) return(degenerate)
  if (decomposition$rank < ncol(design)) {
    kept <- sort(decomposition$pivot[seq_len(decomposition$rank)])
    design <- design[, kept, drop = FALSE]
  }
  
  # Retrieve the coefficient BY NAME so pivoting cannot mis-index it.
  target <- match(coef_name, colnames(design))
  if (is.na(target)) return(degenerate)
  if (nrow(design) <= ncol(design)) return(degenerate)
  
  beta      <- qr.solve(design, y)
  residuals <- as.numeric(y - design %*% beta)
  bread <- tryCatch(solve(crossprod(design)), error = function(e) NULL)
  if (is.null(bread)) return(degenerate)
  
  meat <- matrix(0, ncol(design), ncol(design))
  for (group in unique(cluster)) {
    rows <- cluster == group
    score <- crossprod(design[rows, , drop = FALSE], residuals[rows])
    meat <- meat + score %*% t(score)
  }
  standard_errors <- sqrt(diag(bread %*% meat %*% bread))
  
  list(est  = as.numeric(beta)[target],
       se   = as.numeric(standard_errors)[target],
       rank = ncol(design),
       note = "ok")
}


# ---- TEST BLOCK 7 -> validates SECTION 7 ------------------------------------
local({
  set.seed(RNG_SEED)
  n <- 200L
  x <- rnorm(n)
  cluster_id <- rep(seq_len(20L), each = 10L)
  y <- 2.5 * x + rnorm(n, sd = 0.2)
  design <- cbind(`(int)` = 1, date = x)
  
  fit <- ols_cluster(design, y, cluster_id)
  check_that(identical(fit$note, "ok"), "ols_cluster failed on well-posed data")
  check_that(abs(fit$est - 2.5) < 0.1,
             "ols_cluster did not recover a known slope")
  check_that(is.finite(fit$se) && fit$se > 0,
             "ols_cluster returned an invalid standard error")
  
  # Rank-deficient input must degrade, never error or mis-index.
  collinear <- cbind(design, copy = x)
  check_that(identical(ols_cluster(collinear, y, cluster_id)$note, "ok"),
             "ols_cluster mishandled a collinear column")
  constant <- cbind(`(int)` = 1, date = rep(1, n))
  check_that(identical(ols_cluster(constant, y, cluster_id)$note, "degenerate"),
             "ols_cluster did not flag a constant predictor")
  log_step("TEST 7 passed - cluster-robust estimator recovers a known slope")
})


# =============================================================================
# SECTION 8 - THE DIAGNOSTIC BATTERY
# -----------------------------------------------------------------------------
# Purpose : Estimate the temporal slope for every balance under corrected and
#           uncorrected compositions, adjusted and restricted specifications.
# Inputs  : `P_observed`, `P_ridge`, `P_nnls`, `analysis_data`.
# Outputs : `battery_main`, `battery_nnls`, `battery_within` result tables.
# =============================================================================

# Each coefficient must pass BOTH the full-sample fit and the zero-exclusion
# refit. Detection-limit sensitivity alone is insufficient: the retracted
# woodland decline shrank only from -1.17 to -0.76 across a 75-fold change in
# the detection limit because moving the limit shifts all extremes together.
# Excluding zero-replaced observations is the test that actually discriminates.
build_design <- function(data, adjust) {
  design <- cbind(`(int)` = 1, date = data$date_sc)
  if (adjust %in% c("hab", "habsch")) {
    design <- cbind(design, model.matrix(~ habitat, data)[, -1, drop = FALSE])
  }
  if (adjust %in% c("school", "habsch")) {
    design <- cbind(design, model.matrix(~ school, data)[, -1, drop = FALSE])
  }
  design
}

check_balance <- function(P, totals, which_balance, data,
                          subset = NULL,
                          adjust = c("none", "hab", "habsch", "school"),
                          label = "") {
  adjust <- match.arg(adjust)
  y <- make_balances(P, totals)[[which_balance]]
  
  keep <- if (is.null(subset)) rep(TRUE, nrow(P)) else subset
  keep[is.na(keep)] <- FALSE
  
  design <- build_design(data, adjust)
  full <- ols_cluster(design[keep, , drop = FALSE], y[keep], data$painter[keep])
  
  # Refit excluding artworks whose balance rests on a replaced zero.
  no_zero <- rowSums(P[, BALANCE_PARTS[[which_balance]], drop = FALSE] <= 1e-9) == 0
  keep_nz <- keep & no_zero
  restricted <- if (sum(keep_nz) > MIN_SUBSET_N) {
    ols_cluster(design[keep_nz, , drop = FALSE], y[keep_nz],
                data$painter[keep_nz])
  } else {
    list(est = NA_real_, se = NA_real_, note = "n<min")
  }
  
  spearman <- suppressWarnings(
    stats::cor(y[keep], data$date_sc[keep], method = "spearman",
               use = "complete.obs")
  )
  survives <- is.finite(full$est) && is.finite(restricted$est) &&
    is.finite(restricted$se) && full$est * restricted$est > 0 &&
    (restricted$est - 1.96 * restricted$se) *
    (restricted$est + 1.96 * restricted$se) > 0
  
  data.frame(
    label     = label,
    balance   = which_balance,
    n         = as.integer(sum(keep)),
    est       = as.numeric(full$est)[1],
    lo        = as.numeric(full$est - 1.96 * full$se)[1],
    hi        = as.numeric(full$est + 1.96 * full$se)[1],
    n_nonzero = as.integer(sum(keep_nz)),
    est_nz    = as.numeric(restricted$est)[1],
    lo_nz     = as.numeric(restricted$est - 1.96 * restricted$se)[1],
    hi_nz     = as.numeric(restricted$est + 1.96 * restricted$se)[1],
    spearman  = as.numeric(spearman)[1],
    survives  = isTRUE(survives),
    note      = paste(full$note, restricted$note, sep = "/"),
    stringsAsFactors = FALSE
  )
}

dutch_core_subset <- analysis_data$dutch_core

BATTERY_SPECS <- list(
  list(P = "P_observed", subset = NULL,                adjust = "none",   label = "A uncorrected, unadjusted"),
  list(P = "P_observed", subset = NULL,                adjust = "hab",    label = "B uncorrected, + expert scene type"),
  list(P = "P_observed", subset = NULL,                adjust = "habsch", label = "B uncorrected, + scene type + school"),
  list(P = "P_observed", subset = "dutch_core_subset", adjust = "none",   label = "C uncorrected, Dutch core"),
  list(P = "P_ridge",    subset = NULL,                adjust = "none",   label = "A corrected, unadjusted"),
  list(P = "P_ridge",    subset = NULL,                adjust = "hab",    label = "B corrected, + expert scene type"),
  list(P = "P_ridge",    subset = NULL,                adjust = "habsch", label = "B corrected, + scene type + school"),
  list(P = "P_ridge",    subset = "dutch_core_subset", adjust = "none",   label = "C corrected, Dutch core")
)

battery_main <- bind_rows(lapply(RESPONSES, function(b) {
  bind_rows(lapply(BATTERY_SPECS, function(spec) {
    check_balance(get(spec$P), part_totals, b, analysis_data,
                  subset = if (is.null(spec$subset)) NULL else get(spec$subset),
                  adjust = spec$adjust, label = spec$label)
  }))
}))

# The failure mode kept for the discussion: NNLS manufactures a "credible"
# effect that collapses once zero-replaced artworks are excluded.
battery_nnls <- bind_rows(
  check_balance(P_nnls, part_totals, "b3", analysis_data, NULL, "habsch",
                "NNLS corrected, adjusted"),
  check_balance(P_nnls, part_totals, "b3", analysis_data, dutch_core_subset,
                "none", "NNLS corrected, Dutch core")
)

# Within externally defined scene types: does composition move inside a fixed
# scene? Apparently credible effects appear and collapse under zero exclusion.
battery_within <- bind_rows(lapply(levels(analysis_data$habitat), function(h) {
  in_scene <- analysis_data$habitat == h
  if (sum(in_scene) < MIN_SUBSET_N) return(NULL)
  bind_rows(lapply(RESPONSES, function(b) {
    check_balance(P_observed, part_totals, b, analysis_data, in_scene, "none",
                  sprintf("within %s (n=%d)", h, sum(in_scene)))
  }))
}))

log_step(sprintf("diagnostic battery: %d rows, %d flagged as surviving",
                 nrow(battery_main), sum(battery_main$survives)))


# ---- TEST BLOCK 8 -> validates SECTION 8 ------------------------------------
local({
  expected_rows <- length(RESPONSES) * length(BATTERY_SPECS)
  check_that(nrow(battery_main) == expected_rows,
             sprintf("battery has %d rows, expected %d",
                     nrow(battery_main), expected_rows))
  check_that(nrow(battery_nnls) == 2L, "NNLS battery has unexpected size")
  check_that(nrow(battery_within) > 0L, "within-scene battery produced no rows")
  
  check_that(all(battery_main$n > 0L), "battery rows with empty subsets")
  check_that(all(battery_main$n <= EXPECTED_N), "battery subset larger than corpus")
  check_that(all(battery_main$n_nonzero <= battery_main$n),
             "zero-excluded subset larger than its parent")
  check_that(!anyNA(battery_main$est),
             "diagnostic battery returned NA point estimates")
  check_that(all(battery_main$lo <= battery_main$hi),
             "confidence bounds are inverted")
  check_that(all(battery_main$spearman >= -1 & battery_main$spearman <= 1),
             "Spearman correlation outside [-1, 1]")
  check_that(is.logical(battery_main$survives), "survives must be logical")
  
  dutch_rows <- battery_main[grepl("Dutch core", battery_main$label), ]
  check_that(all(dutch_rows$n == EXPECTED_DUTCH_CORE),
             "restricted rows do not use the Dutch core subset")
  log_step("TEST 8 passed - diagnostic battery internally consistent")
})


# =============================================================================
# SECTION 9 - VARIANCE PARTITION
# -----------------------------------------------------------------------------
# Purpose : Compare the share of compositional variance explained by year,
#           school, painter identity and expert scene type.
# Inputs  : `balances_observed`, `analysis_data`.
# Outputs : `variance_partition` (4 predictors x 4 balances matrix of R-squared).
# =============================================================================

# NOTE: these are single-predictor fixed-effect R-squared values and are
# INFLATED for painter identity, which contributes ~100 dummy columns. They are
# a descriptive screen only. The reportable quantities are the multilevel
# intraclass correlations from the Stage 2 models (section 11), which are
# substantially smaller. Do not quote these numbers as variance components.
VARIANCE_PREDICTORS <- c(year = "date_sc", school = "school",
                         painter = "painter", scene_type = "habitat")

variance_partition <- vapply(balances_observed, function(y) {
  vapply(VARIANCE_PREDICTORS, function(predictor) {
    summary(stats::lm(y ~ analysis_data[[predictor]]))$r.squared
  }, numeric(1))
}, numeric(length(VARIANCE_PREDICTORS)))

school_by_century <- table(analysis_data$school,
                           as_century(analysis_data$median_date))

works_per_painter <- analysis_data |>
  count(painter, name = "n_works") |>
  arrange(desc(n_works))

print(round(variance_partition, 3))
print(school_by_century)


# ---- TEST BLOCK 9 -> validates SECTION 9 ------------------------------------
local({
  check_that(identical(dim(variance_partition),
                       c(length(VARIANCE_PREDICTORS), length(RESPONSES))),
             "variance partition has unexpected dimensions")
  expect_finite(variance_partition)
  check_that(all(variance_partition >= 0 & variance_partition <= 1),
             "R-squared outside [0, 1]")
  check_that(all(variance_partition["painter", ] >= variance_partition["year", ]),
             "painter explains less than year; check the model matrices")
  check_that(sum(school_by_century) == EXPECTED_N,
             "school x century table does not sum to the corpus size")
  check_that(sum(works_per_painter$n_works) == EXPECTED_N,
             "works-per-painter does not sum to the corpus size")
  log_step("TEST 9 passed - variance partition and grouping tables consistent")
})


# =============================================================================
# SECTION 10 - PAINTING-LEVEL BOOTSTRAP OF THE CORRECTION MATRIX
# -----------------------------------------------------------------------------
# Purpose : Propagate uncertainty in the confusion matrix, which is estimated
#           from a small number of validation artworks, into the temporal slopes.
# Inputs  : per_painting_confusion_long.csv, `P_observed`, `analysis_data`.
# Outputs : `bootstrap_draws` (N_BOOTSTRAP x 4) and `bootstrap_quantiles`.
# =============================================================================

# Pixel-level noise is negligible at millions of pixels per artwork; the
# uncertainty that matters is which artworks happened to be annotated.
per_painting_confusion <- read_csv(resolve_input(FILE_PER_PAINT),
                                   show_col_types = FALSE)
expect_cols(per_painting_confusion, c("stem", "true_class", "pred_class", "n"))

class_to_part <- unlist(lapply(PART_NAMES, function(part) {
  setNames(rep(part, length(PARTS[[part]])), PARTS[[part]])
}))

confusion_pooled <- per_painting_confusion |>
  filter(true_class %in% names(class_to_part),
         pred_class %in% names(class_to_part)) |>
  mutate(true_part = class_to_part[true_class],
         pred_part = class_to_part[pred_class]) |>
  group_by(stem, true_part, pred_part) |>
  summarise(n = sum(n), .groups = "drop")

validation_stems <- unique(confusion_pooled$stem)
log_step(sprintf("bootstrap pool: %d validation artworks",
                 length(validation_stems)))

#' Rebuild a row-normalised 5x5 confusion matrix from a sample of artworks.
confusion_from_stems <- function(selected_stems) {
  M <- matrix(0, 5, 5, dimnames = list(PART_NAMES, PART_NAMES))
  for (stem_id in selected_stems) {
    rows <- confusion_pooled[confusion_pooled$stem == stem_id, ]
    M[cbind(rows$true_part, rows$pred_part)] <-
      M[cbind(rows$true_part, rows$pred_part)] + rows$n
  }
  M / pmax(rowSums(M), 1)
}

bootstrap_design <- build_design(analysis_data, "habsch")

bootstrap_one_draw <- function(draw_index) {
  tryCatch({
    resampled <- confusion_from_stems(
      sample(validation_stems, length(validation_stems), replace = TRUE)
    )
    balances <- make_balances(
      correct_composition(P_observed, resampled, "ridge"), part_totals
    )
    vapply(RESPONSES, function(b) {
      as.numeric(ols_cluster(bootstrap_design, balances[[b]],
                             analysis_data$painter)$est)[1]
    }, numeric(1))
  }, error = function(e) {
    warning(sprintf("bootstrap draw %d failed: %s", draw_index,
                    conditionMessage(e)), call. = FALSE)
    setNames(rep(NA_real_, length(RESPONSES)), RESPONSES)
  })
}

set.seed(RNG_SEED)
bootstrap_draws <- t(vapply(seq_len(N_BOOTSTRAP), bootstrap_one_draw,
                            numeric(length(RESPONSES))))
colnames(bootstrap_draws) <- RESPONSES

bootstrap_quantiles <- t(apply(bootstrap_draws, 2, stats::quantile,
                               probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
n_usable_draws <- sum(stats::complete.cases(bootstrap_draws))

log_step(sprintf("bootstrap: %d of %d draws usable", n_usable_draws,
                 N_BOOTSTRAP))
print(round(bootstrap_quantiles, 3))


# ---- TEST BLOCK 10 -> validates SECTION 10 ----------------------------------
local({
  check_that(length(validation_stems) > 1L,
             "bootstrap needs more than one validation artwork")
  check_that(identical(dim(bootstrap_draws),
                       c(N_BOOTSTRAP, length(RESPONSES))),
             "bootstrap draw matrix has unexpected dimensions")
  check_that(n_usable_draws >= 0.9 * N_BOOTSTRAP,
             sprintf("only %d of %d bootstrap draws succeeded",
                     n_usable_draws, N_BOOTSTRAP))
  expect_finite(bootstrap_quantiles)
  check_that(all(bootstrap_quantiles[, 1] <= bootstrap_quantiles[, 3]),
             "bootstrap interval bounds are inverted")
  
  reference <- confusion_from_stems(validation_stems)
  check_that(max(abs(rowSums(reference) - 1)) < 1e-9,
             "rebuilt confusion matrix rows do not normalise")
  log_step("TEST 10 passed - bootstrap complete and well formed")
})


# =============================================================================
# SECTION 11 - STAGE 2: BAYESIAN MULTIVARIATE MODELS
# -----------------------------------------------------------------------------
# Purpose : Fit the four balances jointly with painter (and school) random
#           effects, plus spline and Student-t variants, and compare by PSIS-LOO.
# Inputs  : `analysis_data`, `P_observed` or `P_ridge` via COMPOSITION.
# Outputs : `model_fits`, `loo_results`, `loo_table`, `posterior_summaries`.
# =============================================================================

# COMPOSITION selects which composition feeds the models:
#   "obs"   uncorrected   -> PRIMARY. What the pipeline actually measures.
#   "ridge" bias-corrected -> SENSITIVITY. The correction is estimated from a
#           small validation set and materially reorders the balances
#           (r = 0.65-0.84 against uncorrected), so it cannot carry the headline.
RUN_STAGE2 <- TRUE
COMPOSITION <- Sys.getenv("PAINTINGS_COMPOSITION", "obs")
check_that(COMPOSITION %in% c("obs", "ridge"),
           "PAINTINGS_COMPOSITION must be 'obs' or 'ridge'")

P_for_models <- switch(COMPOSITION, obs = P_observed, ridge = P_ridge)

# Attach the balances, replacing any earlier copies rather than suffixing them.
model_data <- analysis_data[, !grepl("^b[1-4]$", names(analysis_data)),
                            drop = FALSE]
model_data <- bind_cols(model_data, make_balances(P_for_models, part_totals))

log_step(sprintf("balances attached from '%s' composition; n = %d",
                 COMPOSITION, nrow(model_data)))

model_fits          <- list()
loo_results         <- list()
loo_table           <- NULL
posterior_summaries <- NULL
design_analysis     <- NULL

if (RUN_STAGE2) {
  suppressPackageStartupMessages({
    library(brms)
    library(loo)
    library(posterior)
    library(ggplot2)
  })
  options(mc.cores = min(4L, parallel::detectCores()))
  
  # Multivariate brms models require `resp` on every response-specific prior;
  # LKJ priors on rescor and cor are global and must NOT carry `resp`.
  model_priors <- c(
    do.call(c, lapply(RESPONSES, function(response) c(
      set_prior("normal(0, 1)",         class = "b",         resp = response),
      set_prior("student_t(3, 0, 2.5)", class = "Intercept", resp = response),
      set_prior("exponential(1)",       class = "sd",        resp = response),
      set_prior("exponential(1)",       class = "sigma",     resp = response)
    ))),
    set_prior("lkj(2)", class = "rescor"),
    set_prior("lkj(2)", class = "cor")
  )
  spline_priors <- c(
    model_priors,
    do.call(c, lapply(RESPONSES, function(response) {
      set_prior("exponential(1)", class = "sds", resp = response)
    }))
  )
  
  formula_painter <- bf(mvbind(b1, b2, b3, b4) ~ date_sc + lat_sc + lon_sc +
                          habitat + (1 | p | painter))
  formula_school  <- bf(mvbind(b1, b2, b3, b4) ~ date_sc + lat_sc + lon_sc +
                          habitat + (1 | p | painter) + (1 | q | school))
  formula_spline  <- bf(mvbind(b1, b2, b3, b4) ~ s(date_sc) + lat_sc + lon_sc +
                          habitat + (1 | p | painter))
  
  fit_brms <- function(formula, data, prior, adapt_delta = 0.95,
                       family_suffix = NULL) {
    spec <- if (is.null(family_suffix)) {
      formula + set_rescor(TRUE)
    } else {
      formula + family_suffix + set_rescor(TRUE)
    }
    brm(spec, data = data, prior = prior, chains = 4L, iter = 4000L,
        warmup = 2000L, seed = RNG_SEED,
        control = list(adapt_delta = adapt_delta), refresh = 0L)
  }
  
  log_step("fitting M1 (painter random effects)")
  model_fits$M1 <- fit_brms(formula_painter, model_data, model_priors)
  
  log_step("fitting M2 (+ school random effects)")
  model_fits$M2 <- fit_brms(formula_school, model_data, model_priors)
  
  log_step("fitting M3 (restricted to the Dutch core)")
  model_fits$M3 <- fit_brms(formula_painter, filter(model_data, dutch_core),
                            model_priors)
  
  log_step("fitting M4 (thin-plate spline on date)")
  model_fits$M4 <- fit_brms(formula_spline, model_data, spline_priors,
                            adapt_delta = 0.99)
  
  log_step("fitting M5 (Student-t likelihood)")
  model_fits$M5 <- fit_brms(formula_painter, model_data, model_priors,
                            family_suffix = student())
  
  # PSIS-LOO. M3 is fitted to a different number of observations, so it is a
  # restriction to report alongside the others, never a competitor in the table.
  loo_comparable <- c("M1", "M2", "M4", "M5")
  loo_results <- lapply(model_fits[loo_comparable], loo)
  loo_table <- loo_compare(loo_results)
  print(loo_table)
  
  ## --- Posterior summaries -----------------------------------------------
  summarise_fit <- function(fit, model_name) {
    draws <- summarise_draws(
      as_draws_df(fit), mean, sd,
      ~quantile(.x, probs = c(0.025, 0.975)), rhat, ess_bulk, ess_tail
    )
    draws <- draws[!grepl("^r_|^z_|^L_|^lp__|^lprior", draws$variable), ]
    draws$model <- model_name
    draws
  }
  posterior_summaries <- bind_rows(
    lapply(names(model_fits), function(nm) summarise_fit(model_fits[[nm]], nm))
  )
  
  ## --- Posterior predictive checks ---------------------------------------
  for (response in RESPONSES) {
    ggsave(file.path(FIG_DIR, sprintf("ppc_posterior_%s.png", response)),
           pp_check(model_fits$M1, resp = response, ndraws = 100L),
           width = 5, height = 3.5, dpi = 200)
  }
  if (RUN_PRIOR_PREDICTIVE) {
    prior_only <- update(model_fits$M1, sample_prior = "only",
                         seed = RNG_SEED, refresh = 0L)
    ggsave(file.path(FIG_DIR, "ppc_prior_b3.png"),
           pp_check(prior_only, resp = "b3", ndraws = 100L),
           width = 5, height = 3.5, dpi = 200)
  }
  
  ## --- Design analysis ----------------------------------------------------
  # b3 = sqrt(1/2) * log(wooded / open), so halving the wooded:open ratio across
  # the corpus span corresponds to the slope below. Expressing power against an
  # interpretable effect turns a null result into a bound.
  rope_b3 <- sqrt(0.5) * log(0.5) / diff(range(model_data$date_sc))
  log_step(sprintf("halving of wooded:open across the corpus = slope %.3f",
                   rope_b3))
  
  estimate_power <- function(effects = POWER_EFFECTS, nsim = POWER_NSIM) {
    null_model <- lmer(b3 ~ date_sc + (1 | painter), data = model_data)
    sd_painter <- as.data.frame(VarCorr(null_model))$sdcor[1]
    sd_resid   <- sigma(null_model)
    template   <- data.frame(date_sc = model_data$date_sc,
                             painter = model_data$painter)
    
    vapply(effects, function(effect) {
      mean(replicate(nsim, {
        template$y <- effect * template$date_sc +
          rnorm(nlevels(template$painter), 0, sd_painter)[as.integer(template$painter)] +
          rnorm(nrow(template), 0, sd_resid)
        interval <- confint(
          suppressMessages(lmer(y ~ date_sc + (1 | painter), data = template)),
          parm = "date_sc", method = "Wald"
        )
        interval[1] * interval[2] > 0
      }))
    }, numeric(1))
  }
  
  set.seed(RNG_SEED)
  design_analysis <- data.frame(effect = POWER_EFFECTS,
                                power  = estimate_power())
  design_analysis$rope_b3 <- rope_b3
  print(design_analysis)
} else {
  log_step("STAGE 2 skipped (set PAINTINGS_STAGE2=TRUE to fit the brms models)")
}


# ---- TEST BLOCK 11 -> validates SECTION 11 ----------------------------------
local({
  check_that(all(RESPONSES %in% names(model_data)),
             "balances were not attached to the model data")
  expect_finite(model_data[RESPONSES])
  check_that(nrow(model_data) == EXPECTED_N,
             "model data lost rows when balances were attached")
  check_that(sum(model_data$dutch_core) == EXPECTED_DUTCH_CORE,
             "restricted subset drifted before fitting")
  
  if (RUN_STAGE2) {
    check_that(length(model_fits) == 5L, "expected five fitted models")
    check_that(all(vapply(model_fits, inherits, logical(1), "brmsfit")),
               "one or more fits is not a brmsfit object")
    check_that(!is.null(posterior_summaries) && nrow(posterior_summaries) > 0L,
               "posterior summaries are empty")
    expect_cols(posterior_summaries,
                c("variable", "mean", "rhat", "ess_bulk", "ess_tail", "model"))
    
    worst_rhat <- max(posterior_summaries$rhat, na.rm = TRUE)
    min_ess    <- min(posterior_summaries$ess_bulk, na.rm = TRUE)
    check_that(worst_rhat < 1.05,
               sprintf("convergence failure: max R-hat = %.4f", worst_rhat))
    check_that(min_ess > 400,
               sprintf("insufficient effective sample size: min ESS = %.0f",
                       min_ess))
    check_that(!is.null(design_analysis) && all(design_analysis$power >= 0 &
                                                  design_analysis$power <= 1),
               "power estimates outside [0, 1]")
    log_step(sprintf("TEST 11 passed - 5 models converged (max R-hat %.4f)",
                     worst_rhat))
  } else {
    log_step("TEST 11 passed - Stage 1 model data valid (Stage 2 skipped)")
  }
})


# =============================================================================
# SECTION 12 (N-2) - DIAGNOSTICS
# -----------------------------------------------------------------------------
# Purpose : Report convergence, memory, environment and run-flag state so a
#           reviewer can judge the run without re-executing it.
# Inputs  : `posterior_summaries`, session state.
# Outputs : Console diagnostics plus session_info.txt and convergence_summary.csv.
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("DIAGNOSTICS\n")
cat(strrep("=", 78), "\n", sep = "")

## --- Execution flags ---------------------------------------------------------
run_flags <- data.frame(
  flag  = c("RUN_STAGE2", "RUN_PIXEL_COMPARISON", "RUN_PRIOR_PREDICTIVE",
            "COMPOSITION", "RNG_SEED", "N_BOOTSTRAP", "POWER_NSIM",
            "DATA_DIR", "OUT_DIR"),
  value = c(as.character(RUN_STAGE2), as.character(RUN_PIXEL_COMPARISON),
            as.character(RUN_PRIOR_PREDICTIVE), COMPOSITION,
            as.character(RNG_SEED), as.character(N_BOOTSTRAP),
            as.character(POWER_NSIM), DATA_DIR, OUT_DIR),
  stringsAsFactors = FALSE
)
print(run_flags, row.names = FALSE)

## --- Convergence -------------------------------------------------------------
convergence_summary <- NULL
if (RUN_STAGE2 && !is.null(posterior_summaries)) {
  convergence_summary <- posterior_summaries |>
    group_by(model) |>
    summarise(
      n_parameters = n(),
      max_rhat     = max(rhat, na.rm = TRUE),
      min_ess_bulk = min(ess_bulk, na.rm = TRUE),
      min_ess_tail = min(ess_tail, na.rm = TRUE),
      .groups = "drop"
    )
  cat("\nConvergence by model:\n")
  print(as.data.frame(convergence_summary), row.names = FALSE)
  write_csv(convergence_summary, file.path(OUT_DIR, "convergence_summary.csv"))
} else {
  cat("\nConvergence: not applicable (Stage 2 not run).\n")
}

## --- Memory ------------------------------------------------------------------
# gc() returns Ncells/Vcells rows; columns 2 and 6 are the "used" and
# "max used" figures already expressed in MB.
gc_report <- gc(verbose = FALSE)
cat(sprintf("\nMemory in use: %.1f MB (peak %.1f MB)\n",
            sum(gc_report[, 2]), sum(gc_report[, 6])))

global_objects <- ls(envir = globalenv())
if (length(global_objects) > 0L) {
  object_sizes <- vapply(global_objects, function(nm) {
    as.numeric(utils::object.size(get(nm, envir = globalenv())))
  }, numeric(1))
  cat("Largest objects (MB):\n")
  print(round(head(sort(object_sizes, decreasing = TRUE), 5L) / 1024^2, 2))
}

## --- Environment -------------------------------------------------------------
cat(sprintf("\nR version: %s\nPlatform : %s\n",
            R.version.string, R.version$platform))
loaded_versions <- vapply(
  c(REQUIRED_PACKAGES, if (RUN_STAGE2) STAGE2_PACKAGES),
  function(pkg) as.character(utils::packageVersion(pkg)), character(1)
)
print(data.frame(package = names(loaded_versions),
                 version = unname(loaded_versions)), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))

cat(sprintf("\nElapsed: %.1f s\n",
            as.numeric(difftime(Sys.time(), RUN_STARTED_AT, units = "secs"))))


# ---- TEST BLOCK 12 -> validates SECTION 12 ----------------------------------
local({
  check_that(file.exists(file.path(OUT_DIR, "session_info.txt")),
             "session_info.txt was not written")
  check_that(nrow(run_flags) > 0L && !anyNA(run_flags$value),
             "run flags are incomplete")
  if (RUN_STAGE2) {
    check_that(!is.null(convergence_summary) && nrow(convergence_summary) == 5L,
               "convergence summary is missing models")
    check_that(all(convergence_summary$max_rhat < 1.05),
               "a model failed the R-hat threshold")
  }
  log_step("TEST 12 passed - diagnostics captured")
})


# =============================================================================
# SECTION 13 (N-1) - OUTPUT SUMMARY AND ARTIFACT EXPORT
# -----------------------------------------------------------------------------
# Purpose : Format the headline results and write every reportable artifact to
#           OUT_DIR, together with a manifest recording size and row count.
# Inputs  : All result objects produced by sections 2-11.
# Outputs : CSV/TXT artifacts in OUT_DIR and `run_manifest`.
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("OUTPUT SUMMARY\n")
cat(strrep("=", 78), "\n", sep = "")

format_estimate <- function(est, lo, hi, digits = 2) {
  sprintf("%+.*f [%+.*f, %+.*f]", digits, est, digits, lo, digits, hi)
}

HEADLINE_SPEC <- "B uncorrected, + scene type + school"
headline <- battery_main |>
  filter(label == HEADLINE_SPEC) |>
  transmute(balance,
            estimate = format_estimate(est, lo, hi),
            estimate_zero_excluded = format_estimate(est_nz, lo_nz, hi_nz),
            survives)

cat("\nTemporal slope per balance (uncorrected, fully adjusted):\n")
print(as.data.frame(headline), row.names = FALSE)

cat(sprintf("\nBalances surviving the full diagnostic battery: %d of %d\n",
            sum(battery_main$survives), nrow(battery_main)))
cat(sprintf("NNLS demonstration rows flagged as surviving: %d of %d\n",
            sum(battery_nnls$survives), nrow(battery_nnls)))

cat("\nVariance explained (single-predictor R-squared; see section 9 caveat):\n")
print(round(variance_partition, 3))

cat("\nScene-type composition of the corpus:\n")
print(as.data.frame(scene_type_counts), row.names = FALSE)

cat(sprintf("\nGrouping structure: %d artworks, %d painters (%d singletons), %d schools\n",
            nrow(analysis_data), nrow(works_per_painter),
            sum(works_per_painter$n_works == 1L),
            nlevels(analysis_data$school)))

## --- Write artifacts ---------------------------------------------------------
write_csv(battery_main,   file.path(OUT_DIR, "diagnostic_battery.csv"))
write_csv(battery_nnls,   file.path(OUT_DIR, "diagnostic_battery_nnls.csv"))
write_csv(battery_within, file.path(OUT_DIR, "within_scene_type.csv"))
write_csv(works_per_painter, file.path(OUT_DIR, "works_per_painter.csv"))

write_csv(
  as.data.frame(variance_partition) |> tibble::rownames_to_column("predictor"),
  file.path(OUT_DIR, "variance_partition.csv")
)
write_csv(
  as.data.frame(bootstrap_quantiles) |> tibble::rownames_to_column("balance"),
  file.path(OUT_DIR, "bootstrap_quantiles.csv")
)
write_csv(
  as.data.frame.matrix(school_by_century) |> tibble::rownames_to_column("school"),
  file.path(OUT_DIR, "school_by_century.csv")
)

if (RUN_STAGE2) {
  write_csv(posterior_summaries, file.path(OUT_DIR, "posterior_summaries.csv"))
  write_csv(design_analysis,     file.path(OUT_DIR, "design_analysis.csv"))
  capture.output(print(loo_table), file = file.path(OUT_DIR, "loo_comparison.txt"))
}

## --- Manifest ----------------------------------------------------------------
written_files <- unique(c(list.files(OUT_DIR, recursive = TRUE, full.names = FALSE),
                          "run_manifest.csv", "validation_report.csv"))
  run_manifest <- data.frame(
  file       = written_files,
  size_bytes = ifelse(file.exists(file.path(OUT_DIR, written_files)),
                      file.size(file.path(OUT_DIR, written_files)), NA_integer_),
  written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv(run_manifest, file.path(OUT_DIR, "run_manifest.csv"))

cat(sprintf("\nWrote %d artifacts to '%s'\n", nrow(run_manifest), OUT_DIR))


# ---- TEST BLOCK 13 -> validates SECTION 13 ----------------------------------
local({
  check_that(nrow(headline) == length(RESPONSES),
             "headline table does not cover every balance")
  check_that(nrow(run_manifest) > 0L, "manifest is empty")
  check_that(all(run_manifest$size_bytes[!is.na(run_manifest$size_bytes)] > 0),
             "one or more artifacts was written empty")
  log_step("TEST 13 passed - artifacts exported")
})


# =============================================================================
# SECTION 14 (N) - END-TO-END VALIDATION AND TEST SUITE
# -----------------------------------------------------------------------------
# Purpose : Re-read every artifact from disk and re-assert the pipeline's
#           contracts, so a silent failure cannot be mistaken for a clean run.
# Inputs  : Files in OUT_DIR.
# Outputs : Pass/fail report; non-zero exit status on any failure.
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("END-TO-END VALIDATION\n")
cat(strrep("=", 78), "\n", sep = "")

run_test <- function(description, expression) {
  outcome <- tryCatch({
    stopifnot(isTRUE(expression))
    "PASS"
  }, error = function(e) paste0("FAIL: ", conditionMessage(e)))
  data.frame(test = description, result = outcome, stringsAsFactors = FALSE)
}

expected_artifacts <- c(
  "analysis_frame.csv", "school_mapping.csv", "scene_type_counts.csv",
  "diagnostic_battery.csv", "diagnostic_battery_nnls.csv",
  "within_scene_type.csv", "variance_partition.csv",
  "bootstrap_quantiles.csv", "school_by_century.csv",
  "works_per_painter.csv", "session_info.txt", "run_manifest.csv"
)
if (RUN_STAGE2) {
  expected_artifacts <- c(expected_artifacts, "posterior_summaries.csv",
                          "design_analysis.csv", "loo_comparison.txt",
                          "convergence_summary.csv")
}

frame_on_disk   <- read_csv(file.path(OUT_DIR, "analysis_frame.csv"),
                            show_col_types = FALSE)
battery_on_disk <- read_csv(file.path(OUT_DIR, "diagnostic_battery.csv"),
                            show_col_types = FALSE)

test_results <- bind_rows(
  run_test("all expected artifacts exist",
           all(file.exists(file.path(OUT_DIR, expected_artifacts)))),
  run_test("no artifact is zero bytes",
           all(file.size(file.path(OUT_DIR, expected_artifacts)) > 0)),
  run_test("analysis_frame.csv has the expected row count",
           nrow(frame_on_disk) == EXPECTED_N),
  run_test("analysis_frame.csv has no missing values",
           !anyNA(frame_on_disk)),
  run_test("analysis_frame.csv school levels match",
           dplyr::n_distinct(frame_on_disk$school) == EXPECTED_SCHOOL_LEVELS),
  run_test("analysis_frame.csv Dutch core size matches",
           sum(frame_on_disk$dutch_core) == EXPECTED_DUTCH_CORE),
  run_test("analysis_frame.csv artwork identifiers are unique",
           anyDuplicated(frame_on_disk$Image_Name) == 0L),
  run_test("analysis_frame.csv dates fall inside the study period",
           all(frame_on_disk$median_date >= DATE_MIN &
                 frame_on_disk$median_date < DATE_MAX)),
  run_test("diagnostic battery has no missing estimates",
           !anyNA(battery_on_disk$est)),
  run_test("diagnostic battery intervals are ordered",
           all(battery_on_disk$lo <= battery_on_disk$hi)),
  run_test("bootstrap yielded the required proportion of usable draws",
           n_usable_draws >= 0.9 * N_BOOTSTRAP),
  run_test("compositions still close to 1",
           max(abs(rowSums(P_observed) - 1)) < 1e-9),
  run_test("balances are finite",
           all(is.finite(as.matrix(model_data[RESPONSES])))),
  run_test("manifest covers every expected artifact",
           all(expected_artifacts %in% run_manifest$file))
)

if (RUN_STAGE2) {
  posterior_on_disk <- read_csv(file.path(OUT_DIR, "posterior_summaries.csv"),
                                show_col_types = FALSE)
  test_results <- bind_rows(
    test_results,
    run_test("posterior summaries cover all five models",
             dplyr::n_distinct(posterior_on_disk$model) == 5L),
    run_test("all R-hat values are below 1.05",
             max(posterior_on_disk$rhat, na.rm = TRUE) < 1.05),
    run_test("all bulk ESS values exceed 400",
             min(posterior_on_disk$ess_bulk, na.rm = TRUE) > 400),
    run_test("posterior credible intervals are ordered",
             all(posterior_on_disk[["2.5%"]] <= posterior_on_disk[["97.5%"]]))
  )
}

print(test_results, row.names = FALSE)
write_csv(test_results, file.path(OUT_DIR, "validation_report.csv"))

n_failed <- sum(test_results$result != "PASS")
cat("\n", strrep("-", 78), "\n", sep = "")
if (n_failed == 0L) {
  cat(sprintf("PIPELINE COMPLETED SUCCESSFULLY - %d/%d checks passed in %.1f s\n",
              nrow(test_results), nrow(test_results),
              as.numeric(difftime(Sys.time(), RUN_STARTED_AT, units = "secs"))))
} else {
  cat(sprintf("PIPELINE COMPLETED WITH %d FAILING CHECK(S)\n", n_failed))
  print(test_results[test_results$result != "PASS", ], row.names = FALSE)
  if (!interactive()) quit(status = 1L)
  stop("End-to-end validation failed; see validation_report.csv", call. = FALSE)
}
cat(strrep("=", 78), "\n", sep = "")

# =============================================================================
# END OF PIPELINE
# =============================================================================
