# Positive-control Mendelian randomisation for disease-specific mortality GWAS

Code and results for Su C-Y, Butler-Laporte G. *Calibrating the statistical power of disease-specific mortality genome-wide association studies using positive-control Mendelian randomisation* (submitted). An earlier version is on medRxiv: https://doi.org/10.64898/2026.05.29.26354472

## Contents

| Path | What it does |
|---|---|
| `code/1_instruments/` | Format and LD-clump the exposure GWAS (BMI, with the four strength quartiles; household income; within-family BMI) |
| `code/2_mr_mortality/` | Two-sample MR against the six disease-specific mortality GWAS (argument: outcome index 1–6) |
| `code/3_mr_case_control/` | Two-sample MR against the six case-control GWAS (arguments: outcome index, exposure directory, output directory; the same script was run once per exposure) |
| `code/4_figures/build_figures.py` | Figures 1–2 from `results/`: `python code/4_figures/build_figures.py results figures/` |
| `results/{mortality,case_control}/<exposure>/<exposure>--<outcome>.or.txt` | IVW, MR-Egger, weighted median and mode estimates used in the paper |

Paths to local data and software are shown as `<PROJECT_ROOT>`, `<DATA_ROOT>`, `<SOFTWARE_ROOT>` and `<R_LIBS_*>`; set these to your own copies.

## Notes

- In `code/2_mr_mortality/`, a `break` where `next` was intended stops the proxy search, so no proxies were substituted for the mortality outcomes. It is annotated in the code and left unchanged so that the published numbers reproduce.
- Within-family BMI (IEU OpenGWAS `ieu-b-4815`) is on a per-kg/m² scale in `results/`; the figures rescale it to per 1-SD using 4.8 kg/m² (`BMI_SD` in `build_figures.py`).

## Data sources

All inputs are public summary statistics. Disease-specific mortality: https://figshare.com/projects/Progression_GWAS/252002 · BMI: https://portals.broadinstitute.org/collaboration/giant/ · Household income: https://doi.org/10.62891/aac85602 · Within-family BMI: https://opengwas.io/datasets/ieu-b-4815 · Case-control: NIAGADS NG00075, CKDGen round 4, HERMES, GWAS Catalog GCST90129505, GCST90274714 and GCST90104541.

## Software

R 4.1.2 (TwoSampleMR 0.5.6, gwasvcf, data.table, tidyverse), PLINK 1.9, Python ≥ 3.9 (pandas, numpy, matplotlib).
