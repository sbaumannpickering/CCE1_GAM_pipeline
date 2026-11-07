#!/bin/bash
# Reproduce every model, table and figure.
#   CCE1_DATA_DIR / CCE1_RESULTS_DIR override the default data/ and results/ folders.
set -e
cd "$(dirname "$0")/R"

Rscript 02_dvm_models.R              # DVM-amplitude models (six depth bins)
Rscript 01_run_all_models.R          # NASC and cetacean models (DVM as candidate)
Rscript 03_top_day_sensitivity.R     # refit without the 1 and 3 highest days
Rscript 09_cross_validation.R        # leave-one-month-out cross-validation
Rscript 04_build_summary_tables.R    # GAM summary tables (Word/PDF + CSV)
Rscript 05_summary_collages.R        # per-model partial-effect + time-series figures
Rscript 06_dvm_timeseries_figure.R   # DVM amplitude time series by depth bin
Rscript 07_omp_error_figure.R        # OMP closure error
Rscript 08_watermass_regression.R    # water-mass regressions
