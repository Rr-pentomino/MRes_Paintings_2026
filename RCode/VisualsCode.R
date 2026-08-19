#' =============================================================================
#' thesis_figures.R
#' Reproduces every figure in the thesis at publication quality.
#' -----------------------------------------------------------------------------
#' Esmer, A. S. (2026) Computer Vision in Historical Ecology: Quantifying
#' Ecological Baselines Through Automated Semantic Segmentation.
#' MRes Computational Methods in Ecology and Evolution, Imperial College London.
#'
#' Figure 1  Corpus structure in space and time            (map, 3 panels)
#' Figure 2  Painter school against century of production  (heatmap)
#' Figure 3  Segmentation accuracy at two levels           (2 panels)
#' Figure 4  Agreement with human observers                (dot plot)
#' Figure 5  Temporal effect, variance, design sensitivity (3 panels)
#'
#' Every figure is written twice: a vector PDF for submission and a 600 dpi PNG
#' for pasting into Word. A single concatenated All_figures.pdf is also written.
#' No figure carries a title; all explanation belongs in the caption.
#'
#' -----------------------------------------------------------------------------
#' RUN
#'   PAINTINGS_DATA=/path/to/data PAINTINGS_OUT=/path/to/outputs \
#'     Rscript thesis_figures.R
#'   (or set the two paths in the CONFIG block below and press Source in RStudio)
#'
#' INPUTS
#'   $PAINTINGS_OUT/analysis_frame.csv        SECTION 2  of Final STATS.R
#'   $PAINTINGS_OUT/school_by_century.csv     SECTION 13 of Final STATS.R
#'   $PAINTINGS_OUT/posterior_summaries.csv   SECTION 7  of Final STATS.R
#'   $PAINTINGS_OUT/variance_partition.csv    SECTION 9  of Final STATS.R
#'   $PAINTINGS_OUT/design_analysis.csv       SECTION 11 of Final STATS.R
#'   $PAINTINGS_DATA/clean_metadata.xlsx      latitude, longitude
#'
#'   Figures 3 and 4 use per-class IoU and Cohen's kappa values that the Python
#'   notebooks print but do not write to CSV. Those numbers are transcribed in
#'   the FIGURE 3 and FIGURE 4 blocks below, each with its source notebook named,
#'   and are checked against the pooled means reported in the text.
#'
#' PACKAGES
#'   required : readr dplyr tidyr ggplot2 patchwork
#'   Figure 1 : readxl sf rnaturalearth rnaturalearthdata
#'   install.packages(c("readr","dplyr","tidyr","ggplot2","patchwork","readxl",
#'                      "sf","rnaturalearth","rnaturalearthdata"))
#'   If the Figure 1 packages are absent the script reports it and produces the
#'   other four figures rather than stopping.
#' =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(stringr)
})
Sys.setlocale("LC_ALL", "en_GB.UTF-8")
works <- works %>%
  mutate(across(where(is.character), ~ iconv(., from = "", to = "UTF-8", sub = "")))
## =============================================================================
## CONFIG
## =============================================================================
DATA_DIR <- Sys.getenv("PAINTINGS_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("PAINTINGS_OUT",  unset = "outputs")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# cairo_pdf embeds fonts and keeps transparency; fall back if unavailable.
PDF_DEV <- if (isTRUE(unname(capabilities("cairo")))) cairo_pdf else grDevices::pdf

have <- function(...) all(vapply(c(...), requireNamespace, logical(1), quietly = TRUE))

## =============================================================================
## PALETTE AND THEME
## One block, shared by all five figures, so they read as a single system.
## BLUE and ORANGE are distinguishable under deuteranopia, protanopia and
## tritanopia; SEQ_BLUE is single-hue so the heatmap survives greyscale print.
## =============================================================================
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"
SURFACE <- "#fcfcfb"   # panel and page background
INK     <- "#0b0b0b"   # panel tags, darkest text
INK2    <- "#52514e"   # axis titles, category labels
MUTED   <- "#898781"   # tick labels, secondary annotation
GRID    <- "#e1e0d9"   # gridlines, range bars
GREY_MARK <- "#c3c2b7" # unemphasised marks and axis lines
SEA     <- "#eef4fb"; LAND <- "#f0efec"; BORDER <- "#d8d7d0"

SEQ_BLUE <- c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
              "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
              "#0d366b")

#' Minimal theme: no title slot, no panel border, gridlines on one axis only.
theme_thesis <- function(base_size = 9, grid_axis = c("x", "y", "none")) {
  grid_axis <- match.arg(grid_axis)
  th <- theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.35),
      axis.title       = element_text(colour = INK2, size = base_size - 1),
      axis.text        = element_text(colour = MUTED, size = base_size - 1.5),
      axis.line.x      = element_line(colour = GREY_MARK, linewidth = 0.3),
      axis.ticks       = element_blank(),
      legend.title     = element_text(colour = INK2, size = base_size - 1),
      legend.text      = element_text(colour = INK2, size = base_size - 2),
      legend.key       = element_blank(),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      plot.tag         = element_text(colour = INK, size = base_size + 1.5,
                                      face = "bold", hjust = 0),
      plot.tag.position = c(0, 1),
      plot.margin      = margin(6, 8, 4, 4)
    )
  if (grid_axis == "x") th <- th + theme(panel.grid.major.y = element_blank())
  if (grid_axis == "y") th <- th + theme(panel.grid.major.x = element_blank())
  if (grid_axis == "none") th <- th + theme(panel.grid.major = element_blank())
  th
}

#' Place a legend inside the panel. ggplot2 renamed this argument at 3.5.0, so
#' both spellings are handled and the script runs on either side of that break.
legend_inside <- function(x, y, just = c(0, 1)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(x, y),
          legend.justification = just)
  } else {
    theme(legend.position = c(x, y), legend.justification = just)
  }
}

#' Put a colourbar title on the right. Same version split as above.
legend_title_right <- function() {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.title.position = "right")
  } else {
    guides(fill = guide_colourbar(title.position = "right"))
  }
}

FIGS <- list()  # collected for the concatenated PDF

#' Write one figure as vector PDF and 600 dpi PNG.
save_fig <- function(p, stem, width, height) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", device = PDF_DEV)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p,
         width = width, height = height, units = "in", dpi = 600, bg = SURFACE)
  FIGS[[stem]] <<- list(plot = p, width = width, height = height)
  message("wrote ", stem, ".{pdf,png}")
  invisible(p)
}


## =============================================================================
## FIGURE 1 - corpus structure in space and time
## -----------------------------------------------------------------------------
## Marker area encodes the number of works at a location. This matters: the 238
## works occupy only 116 distinct coordinate pairs and one site carries 46 of
## them, so a plot of 238 equal markers collapses to about 116 visible marks and
## hides the concentration the figure exists to show. Panel strips carry painter
## counts as well as work counts, because painters are the effective unit of
## replication (Methods 2.6).
## =============================================================================
fig1 <- NULL
if (!have("readxl", "sf", "rnaturalearth", "rnaturalearthdata")) {
  message("SKIPPING Figure 1: install readxl, sf, rnaturalearth, rnaturalearthdata")
} else {
  suppressPackageStartupMessages({ library(readxl); library(sf); library(rnaturalearth) })
  
  meta <- read_excel(file.path(DATA_DIR, "clean_metadata.xlsx")) |>
    select(Image_Name, latitude, longitude)
  
  works <- read_csv(file.path(OUT_DIR, "analysis_frame.csv"), show_col_types = FALSE) |>
    select(Image_Name, median_date, painter) |>
    left_join(meta, by = "Image_Name") |>
    filter(!is.na(latitude), !is.na(longitude)) |>
    mutate(century = cut(median_date, c(1600, 1700, 1800, 1900),
                         labels = c("1600-1699", "1700-1799", "1800-1899"),
                         right = FALSE))
  
  stopifnot(nrow(works) == 238L, !anyNA(works$century))
  
  panel_lab <- works |>
    group_by(century) |>
    summarise(n = n(), n_painters = n_distinct(painter), .groups = "drop") |>
    mutate(label = sprintf("%s\n%d works, %d painters", century, n, n_painters))
  
  stopifnot(identical(panel_lab$n, c(128L, 27L, 83L)))
  
  sites <- works |>
    count(century, latitude, longitude, name = "n_works") |>
    left_join(select(panel_lab, century, label), by = "century") |>
    mutate(label = factor(label, levels = panel_lab$label))
  
  europe <- ne_countries(scale = "medium", returnclass = "sf")
  XLIM <- c(-12, 28); YLIM <- c(36, 66)   # tight to the data, not a fixed box
  
  fig1 <- ggplot() +
    geom_sf(data = europe, fill = LAND, colour = BORDER, linewidth = 0.18) +
    geom_point(data = sites, aes(longitude, latitude, size = n_works),
               shape = 21, fill = BLUE, colour = "#ffffff",
               stroke = 0.35, alpha = 0.80) +
    facet_wrap(~ label, nrow = 1) +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    scale_size_area(max_size = 7, breaks = c(1, 5, 15, 40),
                    name = "works per location") +
    scale_x_continuous(breaks = seq(-10, 25, 10),
                       labels = function(x) paste0(abs(x), "??", ifelse(x < 0, "W", "E"))) +
    scale_y_continuous(breaks = seq(40, 65, 10),
                       labels = function(y) paste0(y, "??N")) +
    labs(x = NULL, y = NULL) +
    theme_thesis(9, "none") +
    theme(
      panel.background = element_rect(fill = SEA, colour = NA),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.18),
      panel.border     = element_rect(fill = NA, colour = BORDER, linewidth = 0.3),
      panel.spacing    = unit(0.5, "lines"),
      axis.line.x      = element_blank(),
      strip.text       = element_text(colour = INK, size = 8.5, lineheight = 1.15,
                                      margin = margin(b = 4)),
      axis.text        = element_text(colour = MUTED, size = 7),
      legend.position  = "bottom",
      legend.margin    = margin(t = -4)
    )
  
  save_fig(fig1, "Figure1_corpus_map", 6.5, 2.9)
  message(sprintf("  %d works at %d distinct locations; largest site holds %d",
                  nrow(works), nrow(distinct(works, latitude, longitude)),
                  max(sites$n_works)))
}


## =============================================================================
## FIGURE 2 - painter school against century of production
## -----------------------------------------------------------------------------
## School and period are close to collinear in this corpus, which is the
## structural feature limiting identification of H3. Fill is the within-row
## share, so each school's concentration in one period is directly visible;
## cell labels are raw counts; rows are ordered by weighted mean century so the
## near block-diagonal structure reads along the leading diagonal.
## =============================================================================
CENTURY_LABELS <- c("17th", "18th", "19th")

counts <- read_csv(file.path(OUT_DIR, "school_by_century.csv"), show_col_types = FALSE)
stopifnot(all(CENTURY_LABELS %in% names(counts)),
          sum(as.matrix(counts[, CENTURY_LABELS])) == 238L)

long <- counts |>
  pivot_longer(all_of(CENTURY_LABELS), names_to = "century", values_to = "n") |>
  group_by(school) |>
  mutate(total    = sum(n),
         pct      = 100 * n / total,
         mean_cen = sum(n * match(century, CENTURY_LABELS)) / sum(n)) |>
  ungroup() |>
  mutate(school  = factor(school, levels = unique(school[order(-mean_cen, school)])),
         century = factor(century, levels = CENTURY_LABELS,
                          labels = paste0(CENTURY_LABELS, "\ncentury")))

row_totals <- distinct(long, school, total)

fig2 <- ggplot(long, aes(century, school, fill = pct)) +
  geom_tile(colour = SURFACE, linewidth = 1.1) +
  geom_text(aes(label = n, colour = pct > 55,
                fontface = ifelse(n > 0, "bold", "plain")),
            size = 3.1, show.legend = FALSE) +
  geom_text(data = row_totals,
            aes(x = length(CENTURY_LABELS) + 0.72, y = school,
                label = paste("n =", total)),
            inherit.aes = FALSE, hjust = 0, size = 2.9, colour = MUTED) +
  scale_fill_gradientn(colours = SEQ_BLUE, limits = c(0, 100), breaks = c(0, 50, 100),
                       name = "% of school's works") +
  scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = INK)) +
  coord_cartesian(xlim = c(0.5, length(CENTURY_LABELS) + 1.4), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_thesis(9, "none") +
  theme(
    axis.text        = element_text(colour = INK2, size = 8),
    axis.line.x      = element_blank(),
    legend.title      = element_text(colour = INK2, size = 7.5, angle = -90),
    legend.text       = element_text(colour = MUTED, size = 7),
    legend.key.width  = unit(0.22, "cm"),
    legend.key.height = unit(1.9, "cm"),
    plot.margin       = margin(6, 6, 6, 6)
  ) +
  legend_title_right()

save_fig(fig2, "Figure2_school_by_century", 6.3, 3.5)


## =============================================================================
## FIGURE 3 - segmentation accuracy at two levels of aggregation
## -----------------------------------------------------------------------------
## Source: Inference_validation.ipynb, the 16-painting evaluation cohort.
## Panel a is the 15-class scheme with sky and background removed, because
## neither is an ecological reporting class and sky alone (96.5%) would carry
## the mean. Open circles mark classes present in fewer than five works, where
## the IoU is an estimate from very little evidence.
## Panel b is the five-part ecological scheme actually analysed, with 2,000-draw
## bootstrap intervals. The contrast between the two panels is the point of the
## figure: aggregation to the reporting scale nearly doubles accuracy.
## =============================================================================
raw <- tibble::tribble(
  ~cls,             ~iou,  ~n_works,
  "background",     42.65, 16L,
  "building",       60.24, 13L,
  "earth",          25.81, 11L,
  "grass",          32.88, 10L,
  "animal",         63.03, 10L,
  "mountain",       61.95,  7L,
  "path/road",       2.22,  6L,
  "person",         57.68, 13L,
  "rock",            8.97,  8L,
  "shrub/bush",     20.96, 11L,
  "sky",            96.45, 16L,
  "conical tree",    2.63,  3L,
  "broadleaf tree", 49.92, 13L,
  "wooded mass",    14.62,  8L,
  "water",          64.85,  8L
)

analysed <- raw |>
  filter(!cls %in% c("sky", "background")) |>
  arrange(iou) |>
  mutate(cls = factor(cls, levels = cls),
         sparse = n_works < 5)
MEAN_RAW <- mean(analysed$iou)

pooled <- tibble::tribble(
  ~part,              ~iou,  ~lo,  ~hi,
  "wooded",           87.03, 81.0, 96.5,
  "water",            81.03, 51.7, 95.4,
  "human influence",  65.90, 48.2, 82.3,
  "substrate",        59.37, 38.6, 70.8,
  "open vegetation",  39.69, 14.1, 72.9
) |>
  arrange(iou) |>
  mutate(part = factor(part, levels = part))
MEAN_POOLED <- mean(pooled$iou)

# Fail loudly if the transcribed values stop matching the numbers in the text.
stopifnot(abs(MEAN_RAW - 35.83) < 0.05, abs(MEAN_POOLED - 66.60) < 0.05)

f3a <- ggplot(analysed, aes(iou, cls)) +
  geom_segment(aes(x = 0, xend = iou, yend = cls), colour = GREY_MARK, linewidth = 0.6) +
  geom_vline(xintercept = MEAN_RAW, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_point(aes(fill = sparse), shape = 21, size = 2.1, colour = MUTED, stroke = 0.6) +
  annotate("text", x = MEAN_RAW + 2.5, y = 0.9, hjust = 0, size = 2.6, colour = INK2,
           label = sprintf("mean %.1f%%", MEAN_RAW)) +
  scale_fill_manual(values = c(`FALSE` = MUTED, `TRUE` = SURFACE), guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "IoU (%)", y = NULL) +
  theme_thesis(9, "x") +
  theme(axis.text.y = element_text(colour = INK2, size = 7.5))

f3b <- ggplot(pooled, aes(iou, part)) +
  geom_segment(aes(x = lo, xend = hi, yend = part), colour = "#9ec5f4", linewidth = 0.9) +
  geom_vline(xintercept = MEAN_POOLED, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_point(size = 2.6, shape = 21, fill = BLUE, colour = "white", stroke = 0.7) +
  annotate("text", x = MEAN_POOLED - 2.5, y = 0.72, hjust = 1, size = 2.6, colour = INK2,
           label = sprintf("mean %.1f%%", MEAN_POOLED)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "IoU (%)", y = NULL) +
  theme_thesis(9, "x") +
  theme(axis.text.y = element_text(colour = INK2, size = 7.5))

fig3 <- (f3a | f3b) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA)))

save_fig(fig3, "Figure3_segmentation_accuracy", 6.5, 2.85)


## =============================================================================
## FIGURE 4 - agreement with human observers
## -----------------------------------------------------------------------------
## Source: Mask2former.ipynb, quadratic-weighted Cohen's kappa on the
## 20-painting validation cohort. Observer A produced the independent scene
## classification; observer B is the author (Methods 2.5).
## Groups are ordered by the human baseline, so the reader can see at a glance
## where the model matches two humans and where it does not. The grey bar spans
## the three estimates for each group and exists only to lead the eye across.
## Bare soil is the diagnostic case: the observers disagree systematically, so
## no model can agree with both, and soil is excluded from inference.
## =============================================================================
kap <- tibble::tribble(
  ~grp,                 ~hh,    ~mA,    ~mB,
  "water",             0.905,  0.636,  0.675,
  "wooded",            0.690,  0.619,  0.842,
  "total vegetation",  0.641,  0.692,  0.537,
  "animals",           0.600,  0.386,  0.267,
  "human influence",   0.559,  0.481,  0.243,
  "open vegetation",   0.228,  0.275,  0.406,
  "bare soil",        -0.029,  0.626, -0.129
)
POOLED_K <- c(hh = 0.577, mA = 0.684, mB = 0.574)   # quoted in the caption

kap <- kap |>
  mutate(grp = factor(grp, levels = rev(grp)),
         lo  = pmin(hh, mA, mB), hi = pmax(hh, mA, mB))

SERIES <- c(hh = "observer A vs B  (human baseline)",
            mA = "model vs observer A",
            mB = "model vs observer B")

kap_long <- kap |>
  pivot_longer(c(hh, mA, mB), names_to = "series", values_to = "kappa") |>
  mutate(series = factor(SERIES[series], levels = SERIES))

BENCH <- tibble(x = c(0.2, 0.4, 0.6, 0.8),
                lab = c("fair", "moderate", "substantial", "almost perfect"))
NTOP <- nlevels(kap$grp)

fig4 <- ggplot(kap_long, aes(kappa, grp)) +
  geom_segment(data = kap, aes(x = lo, xend = hi, y = grp, yend = grp),
               inherit.aes = FALSE, colour = GRID, linewidth = 1.6) +
  geom_vline(data = BENCH, aes(xintercept = x), colour = GRID, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREY_MARK, linewidth = 0.4) +
  geom_text(data = BENCH, aes(x = x, y = NTOP + 0.55, label = lab),
            inherit.aes = FALSE, size = 2.2, colour = MUTED, vjust = 0) +
  geom_point(aes(shape = series, fill = series, colour = series),
             size = 2.4, stroke = 0.7) +
  scale_shape_manual(values = c(23, 21, 21), name = NULL) +
  scale_fill_manual(values = c(SURFACE, BLUE, ORANGE), name = NULL) +
  scale_colour_manual(values = c(INK2, "white", "white"), name = NULL) +
  scale_x_continuous(limits = c(-0.25, 1.0), breaks = seq(-0.25, 1, 0.25)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  labs(x = "quadratic-weighted Cohen's ??", y = NULL) +
  theme_thesis(9, "x") +
  theme(axis.text.y     = element_text(colour = INK2, size = 8),
        legend.position = "bottom",
        legend.margin   = margin(t = -2),
        plot.margin     = margin(6, 8, 4, 4)) +
  # The same guide must be given to all three aesthetics, otherwise ggplot
  # treats them as different guides and draws the legend two or three times.
  guides(shape  = guide_legend(nrow = 1, override.aes = list(size = 2.4)),
         fill   = guide_legend(nrow = 1, override.aes = list(size = 2.4)),
         colour = guide_legend(nrow = 1, override.aes = list(size = 2.4)))

save_fig(fig4, "Figure4_observer_agreement", 6.4, 3.2)


## =============================================================================
## FIGURE 5 - temporal effect, variance sources, design sensitivity
## -----------------------------------------------------------------------------
## The three panels answer one question in sequence: is there a trend (a), what
## does explain the variation instead (b), and could a trend of the size that
## matters have been detected at all (c).
## Panel y positions are computed by hand rather than dodged, so the offsets are
## explicit and reproducible.
## =============================================================================
post <- read_csv(file.path(OUT_DIR, "posterior_summaries.csv"), show_col_types = FALSE)
BAL    <- c("b1", "b2", "b3", "b4")
BALLAB <- c("b1  natural : human", "b2  vegetation : abiotic",
            "b3  wooded : open", "b4  substrate : water")

#' Posterior mean and 95% credible interval for one variable in one model.
get_post <- function(model, variable) {
  r <- post[post$model == model & post$variable == variable, ]
  stopifnot(nrow(r) == 1L)
  c(mean = r$mean[1], lo = r$`2.5%`[1], hi = r$`97.5%`[1])
}

## --- (a) effect of year -------------------------------------------------------
MODELS <- c(M1 = "M1  scene type", M2 = "M2  + school", M3 = "M3  Dutch only")
OFFSET <- c(M1 = 0.24, M2 = 0.0, M3 = -0.24)

yearfx <- expand.grid(bal = BAL, model = names(MODELS),
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
pull_field <- function(field) {
  mapply(function(m, b) unname(get_post(m, paste0("b_", b, "_date_sc"))[[field]]),
         yearfx$model, yearfx$bal)
}
yearfx$mu <- pull_field("mean")
yearfx$lo <- pull_field("lo")
yearfx$hi <- pull_field("hi")
# y positions are set by hand rather than dodged, so the offsets are explicit:
# b1 sits at 3 and b4 at 0, with the three models offset by +0.24, 0 and -0.24.
yearfx$y  <- (length(BAL) - match(yearfx$bal, BAL)) + unname(OFFSET[yearfx$model])
yearfx$model <- factor(unname(MODELS[yearfx$model]), levels = unname(MODELS))

f5a <- ggplot(yearfx, aes(mu, y, colour = model)) +
  geom_vline(xintercept = 0, colour = INK2, linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.5, alpha = 0.8) +
  geom_point(size = 1.7, shape = 21, aes(fill = model), colour = "white", stroke = 0.5) +
  scale_colour_manual(values = unname(c(BLUE, ORANGE, AQUA)), name = NULL) +
  scale_fill_manual(values = unname(c(BLUE, ORANGE, AQUA)), name = NULL) +
  scale_y_continuous(breaks = rev(seq_along(BAL)) - 1, labels = BALLAB,
                     limits = c(-0.7, length(BAL) + 0.3)) +
  labs(x = "effect of year (per SD, 95 y)", y = NULL) +
  theme_thesis(9, "x") +
  theme(axis.text.y       = element_text(colour = INK2, size = 7.2),
        legend.background = element_blank(),
        legend.text       = element_text(size = 6.6),
        legend.key.height = unit(0.34, "cm")) +
  legend_inside(0.02, 0.99, c(0, 1))

## --- (b) variance sources -----------------------------------------------------
## Intraclass correlation from the group-level SDs: painter from M1, school from
## M2 (the only model containing a school term). Year is the single-predictor
## share from variance_partition.csv, kept on the same axis for comparison only.
icc <- function(model, b, term) {
  sd_term <- get_post(model, sprintf("sd_%s__%s_Intercept", term, b))[["mean"]]
  denom <- get_post(model, sprintf("sd_painter__%s_Intercept", b))[["mean"]]^2 +
    get_post(model, sprintf("sigma_%s", b))[["mean"]]^2
  if (model == "M2") denom <- denom + get_post(model, sprintf("sd_school__%s_Intercept", b))[["mean"]]^2
  100 * sd_term^2 / denom
}

vp <- read_csv(file.path(OUT_DIR, "variance_partition.csv"), show_col_types = FALSE)
vp_year <- vapply(BAL, function(b) 100 * vp[[b]][vp$predictor == "year"], numeric(1))

HGT <- 0.24
varsrc <- bind_rows(
  tibble(bal = BAL, source = "painter", value = vapply(BAL, icc, numeric(1), model = "M1", term = "painter"), dy =  HGT),
  tibble(bal = BAL, source = "school",  value = vapply(BAL, icc, numeric(1), model = "M2", term = "school"),  dy =  0),
  tibble(bal = BAL, source = "year",    value = unname(vp_year),                                              dy = -HGT)
) |>
  mutate(source = factor(source, levels = c("painter", "school", "year")),
         y = (length(BAL) - match(bal, BAL)) + dy)

# The x axis is fixed at 0-30% so the three panels of Figure 5 stay comparable;
# fail loudly rather than silently clipping a bar if a component grows past it.
stopifnot(max(varsrc$value) < 30)

f5b <- ggplot(varsrc, aes(value, y, fill = source)) +
  geom_col(orientation = "y", width = HGT) +
  # Bars below 1% are sub-pixel on a 0-30% axis, so they carry a direct label.
  # Without it a near-zero share is indistinguishable from a missing bar.
  geom_text(data = filter(varsrc, value < 1),
            aes(x = 0.5, y = y, label = sprintf("%.2f%%", value)),
            hjust = 0, size = 2.0, colour = MUTED, inherit.aes = FALSE) +
  scale_fill_manual(values = c(painter = BLUE, school = ORANGE, year = GREY_MARK),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = rev(seq_along(BAL)) - 1, labels = BAL,
                     limits = c(-0.55, length(BAL) - 0.25)) +
  labs(x = "variance explained (%)", y = NULL) +
  theme_thesis(9, "x") +
  theme(axis.text.y       = element_text(colour = INK2, size = 7.5),
        legend.background = element_blank(),
        legend.text       = element_text(size = 6.6),
        legend.key.height = unit(0.34, "cm")) +
  legend_inside(0.99, 0.99, c(1, 1))

## --- (c) design sensitivity ---------------------------------------------------
da   <- read_csv(file.path(OUT_DIR, "design_analysis.csv"), show_col_types = FALSE)
ROPE <- abs(da$rope_b3[1])                       # slope equal to a halving of b3
PW   <- 100 * approx(da$effect, da$power, xout = ROPE)$y

f5c <- ggplot(da, aes(effect, 100 * power)) +
  geom_hline(yintercept = 80, colour = GREY_MARK, linewidth = 0.35, linetype = "22") +
  annotate("text", x = 0.5, y = 82, label = "80%", hjust = 1, vjust = 0,
           size = 2.3, colour = MUTED) +
  geom_line(colour = BLUE, linewidth = 0.7) +
  geom_point(colour = BLUE, fill = BLUE, shape = 21, size = 1.2, stroke = 0.4) +
  # The label sits in the empty upper left, clear of the curve, and the leader
  # starts below the text block so the two never overlap.
  annotate("segment", x = 0.10, xend = ROPE - 0.004, y = 56, yend = PW + 3.5,
           colour = ORANGE, linewidth = 0.3) +
  annotate("point", x = ROPE, y = PW, shape = 21, size = 2.6,
           fill = SURFACE, colour = ORANGE, stroke = 0.8) +
  annotate("text", x = 0.015, y = 74, hjust = 0, vjust = 0.5, size = 2.3,
           colour = ORANGE, lineheight = 1.05,
           label = sprintf("halving of\nwooded : open\n%.0f%% power", PW)) +
  scale_x_continuous(limits = c(0, 0.52), breaks = seq(0, 0.5, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "true slope on b3", y = "power (%)") +
  theme_thesis(9, "y")

fig5 <- (f5a | f5b | f5c) +
  plot_layout(widths = c(1.35, 1.05, 1.0)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA)))

save_fig(fig5, "Figure5_temporal_variance_power", 6.6, 2.95)


## =============================================================================
## CONCATENATED PDF AND CONSOLE REPORT
## =============================================================================
allpdf <- file.path(FIG_DIR, "All_figures.pdf")
PDF_DEV(allpdf, width = 6.6, height = 3.6, onefile = TRUE)
for (f in FIGS) print(f$plot)
invisible(dev.off())
message("wrote All_figures.pdf (", length(FIGS), " figures)")

cat("\n--- values reproduced -----------------------------------------------\n")
cat(sprintf("mean IoU, 13 analysed classes : %.2f%%\n", MEAN_RAW))
cat(sprintf("mean IoU, 5 ecological parts  : %.2f%%\n", MEAN_POOLED))
cat(sprintf("pooled kappa  A / B / A-vs-B  : %.3f / %.3f / %.3f\n",
            POOLED_K[["mA"]], POOLED_K[["mB"]], POOLED_K[["hh"]]))
cat("painter ICC, M1 (%) :", sprintf("%.2f", filter(varsrc, source == "painter")$value), "\n")
cat("school  ICC, M2 (%) :", sprintf("%.2f", filter(varsrc, source == "school")$value), "\n")
cat("year    R2      (%) :", sprintf("%.2f", filter(varsrc, source == "year")$value), "\n")
cat(sprintf("slope for a halving of b3     : %.3f, power %.1f%%\n", ROPE, PW))
cat("\nfigures written to: ", normalizePath(FIG_DIR), "\n", sep = "")
cat("packages attached: ", paste(names(utils::sessionInfo()$otherPkgs), collapse = ", "), "\n", sep = "")
