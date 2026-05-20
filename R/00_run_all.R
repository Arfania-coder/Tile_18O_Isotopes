# =============================================================================
# Master runner: execute the full reproducible analysis pipeline in order.
#
# Usage:
#   * Open Arfania_coder.Rproj (or set the repo root as your working directory)
#   * source("R/00_run_all.R")
#
# This will:
#   1. Run the main analysis (Tables 1-8, Figs 2-4, 6-7, S1-S2)
#   2. Build the Craig diagram with inset (LMWL + evaporation line)
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

scripts <- c(
  here("R", "01_main_analysis.R"),
  here("R", "02_craig_diagram.R")
)

for (s in scripts) {
  cat("\n\n#######################################################\n")
  cat("###  Running:", basename(s), "\n")
  cat("#######################################################\n\n")
  source(s, echo = FALSE)
}

cat("\n\nAll scripts completed. Outputs in:\n")
cat("  ", here("output"),  "\n")
cat("  ", here("figures"), "\n")
