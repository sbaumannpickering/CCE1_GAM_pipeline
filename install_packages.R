# R packages used by the pipeline (versions the published results were produced with).
pkgs <- c("mgcv", "nlme", "R.matlab", "dplyr", "officer", "flextable")
install.packages(setdiff(pkgs, rownames(installed.packages())))
