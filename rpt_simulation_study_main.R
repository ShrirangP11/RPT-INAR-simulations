# Run from the folder containing both R scripts.
# Install glmnet once beforehand: install.packages("glmnet")

source("rpt_simulation_sensitivity.R")

# Global seed: the simulation function does not reset it.
set.seed(12345)

results <- run_rpt_study(
  B = 500,
  N = c(300, 350, 400, 450, 500, 1000),
  Qs = 15
)

results$summary

# Record the R and package versions used for this run.
capture.output(
  sessionInfo(),
  file = file.path("rpt_dft_results", "sessionInfo.txt")
)
