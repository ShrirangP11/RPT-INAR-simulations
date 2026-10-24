# RPT and DFT: Monte Carlo study for seasonal INAR models

R scripts comparing period recovery by the Ramanujan Periodicity Transform
(RPT) and a DFT harmonic score in three seasonal INAR(1) simulation cases.

## Files

- `rpt_simulation_sensitivity.R`: defines the simulation and analysis function
  `run_rpt_study()`. Source this file before calling the function.
- `rpt_simulation_study_main.R`: runs the Monte Carlo comparison with the settings
  below and records the R and package versions.

The simulation function is unchanged from the supplied study code. The launch
script uses relative filenames instead of a computer-specific Downloads path.

## Requirements

R (at least version 3.6.0) and the R package `glmnet`.
No external dataset is needed.

Install the package once, in the R console:

```r
install.packages("glmnet")
```

## Run the main study

Download and extract the repository. Set R's working directory to the extracted
folder containing both scripts. In RStudio, you can open
`rpt_simulation_study_main.R` and select **Session > Set Working Directory >
To Source File Location**.

Then run:

```r
source("rpt_simulation_study_main.R")
```

The launch script sets the global seed to `12345`, as in the supplied launch
file, and uses:

- 500 replications for each case and sample size;
- sample sizes 300, 350, 400, 450, 500, and 1000;
- maximum candidate period Q = 15, with candidates 2 through 15;
- burn-in of 500 observations;
- LASSO was used to estimate the coefficients.
- ten-fold cross-validation, selecting `lambda.min`.

The function uses the existing random-number state. To repeat a run, reset the
global seed before calling it. Keep the full set of settings and software
versions with the final results.

## Simulation cases and recovery criteria

1. **Periodic thinning:** thinning probability
   0.7 * abs(sin(2 * pi * t / 7)), with Poisson innovations of mean 1.
   Recovery requires the highest RPT or DFT score to identify period 7.
2. **Periodic geometric innovations:** constant thinning probability 0.7,
   with geometric success probabilities repeating 0.4, 0.5, 0.6, 0.7, 0.8.
   R's geometric distribution counts failures before the first success.
   Recovery requires the highest score to identify period 5.
3. **Multiple periodic components:** constant thinning probability 0.7,
   with Poisson innovations whose means contain components of periods 3 and 8.
   Recovery requires the two highest scores to identify both 3 and 8.
   The overall combined pattern has period 24.

The known number of components is used for this recovery benchmark.
Both methods analyse the same simulated trajectory in each replication.

RPT scores are sums of absolute fitted coefficients grouped by candidate period.
Zero RPT scores are excluded from selection; an all-zero fit counts as a
failure.

For each candidate period q, the DFT score averages the demeaned periodogram
powers at the nearest Fourier bins to the primitive positive harmonics k/q,
where gcd(k, q) = 1 and k <= floor(q/2). This is a harmonic score, rather than
selection of the single largest periodogram peak. DFT and RPT rank the same
candidate periods 2 through Q.

## Output

The main launch script creates:

- `rpt_dft_results/recovery_percentages.csv`: recovery percentages for both
  methods;
- `rpt_dft_results/sessionInfo.txt`: R and package versions used for the run.

Include the outputs from the final manuscript run in the repository, so readers
can compare their results with the reported table.

## Optional sensitivity analysis for Q

To examine Q = 10, 15, and 20, use:

```r
source("rpt_simulation_sensitivity.R")
set.seed(12345)

sensitivity_results <- run_rpt_study(
  B = 500,
  N = c(300, 350, 400, 450, 500, 1000),
  Qs = c(10, 15, 20),
  outdir = "rpt_dft_Q_sensitivity"
)

capture.output(
  sessionInfo(),
  file = file.path("rpt_dft_Q_sensitivity", "sessionInfo.txt")
)
```

All Q values use the same trajectory and cross-validation folds within each
replication. The main launch script runs Q = 15 only; the optional analysis is
a separate run and is not evidence that Q sensitivity has already been
reported in the manuscript.

## Scope

These files produce Monte Carlo recovery tables. They do not include the
plotting scripts for the original single-trajectory illustrations.
