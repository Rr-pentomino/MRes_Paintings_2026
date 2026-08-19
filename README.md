# Computer Vision in Historical Ecology

Code and derived data for an MRes thesis (Computational Methods in Ecology and
Evolution, Imperial College London) asking whether the land-cover composition of
European landscape paintings changed between 1611 and 1899.

Three semantic segmentation architectures are fine-tuned on manually annotated
oil paintings; the best is applied to 238 works, each of which becomes a
five-part land-cover composition. Because those parts are constrained to sum to
one, the analysis is conducted in isometric log-ratio coordinates rather than on
the proportions themselves, with painter identity as a random effect and an
independently assigned scene type as the only scene covariate.

## Repository

| Folder | Contents |
|---|---|
| `Dataset/Metadata` | `clean_metadata.xlsx` — dates, coordinates, attribution, artistic tradition, and the independent EUNIS-aligned scene classification, `Cohens_kappa_matrix.xlsx` — Ordinal dafor ratings `PaintingsDataset.xlsx` links to National Gallery online archive along with painter and painting names |
| `PythonCode` | Fine-tuning notebooks (`Mask2former`, `DeepLab`, `Segformer`), `Inference_validation` (full-corpus inference and the held-out validation), `Convert_confusion` (bridges Python output to R input) |
| `PythonOutputs` | Per-artwork pixel counts, the pooled confusion matrix, per-artwork confusion matrices in long form, and empirical measurement error |
| `RCode` | `Final STATS.R` — the complete statistical pipeline |
| `ROutputs` | Every table and figure reported in the thesis, plus the run manifest, session info and validation report |

## Pipeline

Run in this order; each step writes what the next one reads.

1. **`PythonCode/Mask2former.ipynb`** (and `DeepLab`, `Segformer` for the
   architecture comparison) — fine-tune on the annotated artworks.
2. **`PythonCode/Inference_validation.ipynb`** — apply the selected model to the
   full corpus and score it against the held-out artworks. Writes
   `composition_counts.csv`, `R_pooled_rownorm.csv`, `measurement_error.csv` and
   `per_painting_confusion.npz`.
3. **`PythonCode/Convert_confusion.ipynb`** — convert the `.npz` to
   `per_painting_confusion_long.csv`. Not optional: R cannot read `.npz`, and the
   painting-level bootstrap in step 4 depends on it.
4. **`RCode/Final STATS.R`** — compositional analysis, diagnostic battery,
   Bayesian models and design analysis.

## Running the analysis

`Final STATS.R` reads all four Python outputs plus `clean_metadata.xlsx` from a
single directory, so copy `clean_metadata.xlsx` alongside the contents of
`PythonOutputs` first. No paths are hardcoded:

```r
Sys.setenv(PAINTINGS_DATA   = "path/to/inputs",
           PAINTINGS_OUT    = "path/to/outputs",
           PAINTINGS_STAGE2 = "TRUE")     # omit for the fast, Stan-free stage
source("RCode/Final STATS.R")
```

Stage 1 (sections 2–10) runs in seconds and reproduces every descriptive and
frequentist result. Stage 2 fits five multivariate `brms` models and takes
considerably longer.

The script validates itself: thirteen inter-section test blocks check dimensions,
missingness and compositional closure as it goes, and a final suite re-reads
every artifact from disk and exits non-zero on failure. `ROutputs/validation_report.csv`
records the outcome of the run that produced the reported numbers, and
`ROutputs/session_info.txt` the exact package versions.

## Data availability

Painting images are reproduced from the National Gallery, London online archive
and are not redistributed here. The derived tables the analysis actually consumes
(pixel counts, curatorial metadata and the scene classification) are included
in full, so every result below the segmentation step is reproducible without
them. Fine-tuned model weights are available on request pending publication.

## Requirements

Python: `torch`, `transformers`, `pillow`, `numpy`, `pandas`, `scipy`.
R: `dplyr`, `readr`, `readxl`, `stringr`, `stringi`, `tidyr`, `tibble`, `nnls`,
`lme4`, `brms`, `loo`, `posterior`, `ggplot2`.

## Licence

Painting images remain under their original terms.

## Author 

Ahmet Selim Esmer // 
esmerahmetselim@gmail.com
