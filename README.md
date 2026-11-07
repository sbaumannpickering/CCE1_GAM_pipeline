# CCE1 automated AR1 GAM pipeline

R code for the statistical analysis (and MATLAB code for the scattering-layer tracing) in Baumann-Pickering et al. (*Limnology and
Oceanography*, in revision): generalized additive models (GAMs) and water-mass regressions with AR1 errors relating active acoustic backscatter (NASC), diel vertical migration (DVM) amplitude and cetacean acoustic occurrence/calling to water-column
oceanography at the California Current Ecosystem mooring CCE1 (Oct 2016 -
Nov 2017).

- **Code archive (DOI):** to be added (Zenodo release of this repository)
- **Data archive (DOI):** [10.5281/zenodo.24044208](https://doi.org/10.5281/zenodo.24044208) (Zenodo)

## Pipeline

For each response (6 NASC depth bins, 6 DVM depth bins, and the cetacean
responses), `00d_gam_pipeline_functions.R` runs the same automated procedure:

1. Every candidate predictor (36 physical variables: temperature, salinity and
   their anomalies at 8 depths, plus dissolved oxygen at 500 m, nitrate and
   chlorophyll a at 40 m, FSLE; cetacean models also get 6 NASC and 6 DVM
   candidates) is fit alone as a smooth `s(x, k = 4, bs = "cs")` with AR1 errors;
   those with p < 0.05 are kept.
2. A pilot fit without correlation gives the lag-1 residual autocorrelation
   used as the AR1 coefficient (`rho`, limited to 0-0.95); `bam(method =
   "fREML", discrete = TRUE)` is used for Gaussian, negative binomial and
   Tweedie responses, `gamm()` (PQL) with `corAR1` for binary occurrence.
3. Candidates are thinned for concurvity (> 0.6), then added by forward
   selection on AIC (minimum improvement 2, at most 6 terms; p-value-based
   for binary models), and pruned after selection (edf < 0.5, p >= 0.05, or
   shrinkage-model p >= 0.05).
4. Final models are checked by refitting without the 1 and 3 highest days
   (`03`) and by leave-one-calendar-month-out cross-validation (`09`; term
   selection is not repeated per fold, so CV scores are optimistic).

DVM amplitude is the night mean minus the day mean of hourly NASC (civil
sunrise/sunset; night = the preceding sunset to that day's sunrise).

## Versions used

MATLAB R2025b with the Statistics and Machine Learning Toolbox is needed only for
the scattering-layer tracing (`matlab/`). R 4.6.0 with mgcv 1.9-4, nlme 3.1-169, R.matlab 3.8.0, dplyr 1.2.1, officer
0.7.7 and flextable 0.10.3 (`install_packages.R`). PDF copies of the Word
tables need Microsoft Word and the `doconv` package (optional; skipped with a
message otherwise). Other mgcv versions can select slightly different terms
when two nearly collinear predictors are close in AIC.

## Running

1. Download the data deposit into `data/` (see `data/README.md`).
2. `./run_all.sh` (may take several hours depending on compute capacity; 
   the binomial `gamm()` fits take longest), or run the scripts in `R/` 
   one at a time from the `R/` directory.
   Set `CCE1_DATA_DIR` / `CCE1_RESULTS_DIR` to use other folders.

Outputs go to `results/` (not tracked by git).

The nightly mesopelagic scattering layer file `MSL_night_seedWalk.csv` is in the
data deposit. To regenerate it from the 5 m WBAT exports, run
`matlab/run_MSL_seedwalk.m` in MATLAB (see `matlab/README.md`); the output goes
to `results/MSL/`. The same folder has the interactive tool for hand-labeling the
layer edge and a script that scores the detector against those labels.

## Scripts, outputs and where they appear in the paper

| Script | Output (in `results/`) | Paper |
|---|---|---|
| `02_dvm_models.R` | `models/DVM*_model.rds`, `dvm_model_summary.csv` | DVM models (Table S7, Fig. S6) |
| `01_run_all_models.R` | `models/<response>_model.rds`, `all_model_summary.csv`, `dvm_term_summary.csv` | NASC and cetacean models (Tables S7, S9; Figs. 6, 7, S7, S9) |
| `03_top_day_sensitivity.R` | `models/top_day_sensitivity.csv` | Table S8 |
| `09_cross_validation.R` | `models/cross_validation.csv` | CV columns of Tables S7, S9 |
| `04_build_summary_tables.R` | `models/table_SA_*.csv`, `table_SB_cetacean_terms.csv`, Word/PDF tables | Tables S7 (S-A), S9 (S-B) |
| `05_summary_collages.R` | `model_plots/<response>_collage.pdf` | Figs. 6, 7, S6, S7, S9 |
| `06_dvm_timeseries_figure.R` | `model_plots/DVM_amplitude_timeseries*.pdf` | Fig. S5 |
| `07_omp_error_figure.R` | `omp_uncertainty/` | Fig. S2, Table S2 |
| `matlab/run_MSL_seedwalk.m` | `MSL/MSL_night_seedWalk.csv` | Input to Fig. 5 and Table S6 (MSL regressions) |
| `matlab/run_MSL_labeling.m`, `matlab/run_MSL_scoring.m` | `MSL/MSL_TopEdge_HandLabels.mat`, `MSL/MSL_score_*.csv` | Hand-labeling tool and detector scoring for the MSL top edge |
| `08_watermass_regression.R` | `fig5_watermass_ar1/` (figure, CSV, Word/PDF table) | Fig. 5, Table S6 (S-C) |

`00a`-`00f` are shared helpers (data loading, DVM calculation, plotting, the
GAM functions, Word-table styling, cetacean model specifications) and
`config.R` holds the paths.

## Not included

- Figures and tables in the paper are assembled from these outputs; some
  (Tables S1, S3, S4, S5, the formatting of S9, Figs. 1-4, S1, S3, S4) are
  not generated by this code.
- The dagger (exploratory) marker in `04` implements only the "fewer than 25
  detection days" criterion. The additional criteria used in the Supporting
  Information (negative cross-validated deviance explained; deviance explained
  falling below half or all terms lost without the 3 highest days) were applied
  by hand from `cross_validation.csv` and `top_day_sensitivity.csv`.
- The MATLAB code that assembled `CCE1_WBAT_HARP_daily.mat` from the CTD,
  nutrient, FSLE and detection records is not part of this release; the file is
  provided in the data deposit. The code that traces the mesopelagic scattering
  layer is in `matlab/`.

## License

Code: MIT (see `LICENSE`).
