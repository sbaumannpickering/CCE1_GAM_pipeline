# Paths shared by every script. All scripts are meant to be run from the R/
# directory (see run_all.sh); defaults are relative to the repository root and
# can be overridden with environment variables.
#
#   CCE1_DATA_DIR     input data (downloaded from the data deposit; see data/README.md)
#   CCE1_RESULTS_DIR  everything the pipeline writes

REPO_ROOT   <- dirname(normalizePath(getwd()))
DATA_DIR    <- Sys.getenv("CCE1_DATA_DIR",    unset = file.path(REPO_ROOT, "data"))
RESULTS_DIR <- Sys.getenv("CCE1_RESULTS_DIR", unset = file.path(REPO_ROOT, "results"))

MODEL_DIR <- file.path(RESULTS_DIR, "models")        # <response>_model.rds, summaries, CSVs
PLOT_DIR  <- file.path(RESULTS_DIR, "model_plots")   # collages and DVM time series
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR,  showWarnings = FALSE, recursive = TRUE)

# input files (names as in the data deposit)
COMBINED_DAILY_MAT <- file.path(DATA_DIR, "combined_diel_daily.mat")
WBAT_EXPORT_MAT    <- file.path(DATA_DIR, "CCE1_SB_01_all.mat")
OMP_DIR            <- file.path(DATA_DIR, "OMP")
MSL_CSV            <- file.path(DATA_DIR, "MSL_final_per_night_seedWalk.csv")
