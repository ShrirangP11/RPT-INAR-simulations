# Simple Monte Carlo comparison: RPT versus a Fourier harmonic period score.
# Set your preferred global random seed in R before running the study.
# The function uses the current random-number state without resetting it.
# Run:
#   source("rpt_simulation_sensitivity.R")
#   results <- run_rpt_study(B = 500)
#   results$summary
#
# Defaults: N=300,350,400,450,500,1000 for EVERY case; burn-in=500; lambda.min only.
# Use Qs=15 if you want just the original candidate bound.
# Recovery: top one period equals 7 or 5 in Cases 1 and 2;
#           top two periods are {3,8} in Case 3.
# The known number of components is used as a simulation benchmark.
#
# DFT_harmonic: for each q, average nearest-bin periodogram powers at
#               primitive positive harmonics k/q (gcd(k,q)=1).
#               Rank candidate periods 2,...,Q, as for RPT.

# Both methods use the SAME trajectory. All Q values use the same CV folds.
# Only a table of recovery percentages is printed and saved.

rpt_gcd <- function(a, b) {
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  a
}

run_rpt_study <- function(B = 500,
    N = c(300, 350, 400, 450, 500, 1000),
    Qs = c(10, 15, 20), burn = 500,
    outdir = "rpt_dft_results") {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop('Install glmnet first: install.packages("glmnet")')
  if (getRversion() < "3.6.0") stop("Use R >= 3.6.")
  stopifnot(length(B) == 1, B >= 1, B == floor(B),
    length(burn) == 1, burn >= 0, burn == floor(burn),
    length(N) >= 1, all(N >= 40),
    all(N == floor(N)), length(Qs) >= 1, all(Qs >= 8),
    all(Qs == floor(Qs)), max(Qs) < min(N))
  N <- sort(unique(N))
  Qs <- sort(unique(Qs))
  max_Q <- max(Qs)
  summary <- data.frame()

  for (case in 1:3) {
    true_periods <- if (case == 1) 7 else if (case == 2) 5 else c(3, 8)
    K <- length(true_periods)
    for (n in N) {
      time <- 1:(n + burn)

      # Build the dictionary once, outside the repeated simulations.
      A <- matrix(numeric(0), nrow = n, ncol = 0)
      column_period <- numeric(0)
      harmonic_bins <- vector("list", max_Q)
      for (q in 2:max_Q) {
        coprimes <- numeric(0)
        for (k in 1:q) {
          if (rpt_gcd(k, q) == 1) coprimes <- c(coprimes, k)
        }
        c_q <- numeric(q)
        for (u in 0:(q - 1)) {
          for (k in coprimes) c_q[u + 1] <- c_q[u + 1] + cos(2 * pi * k * u / q)
        }
        c_q <- round(c_q)
        for (m in 0:(length(coprimes) - 1)) {
          atom <- c_q[((1:n - m) %% q) + 1]
          A <- cbind(A, atom)
          column_period <- c(column_period, q)
        }
        for (k in 1:floor(q / 2)) {
          if (rpt_gcd(k, q) == 1) {
            bin <- floor(n * k / q + 0.5)
            bin <- min(floor(n / 2), max(1, bin))
            harmonic_bins[[q]] <- c(harmonic_bins[[q]], bin)
          }
        }
      }

      rpt_success <- numeric(length(Qs))
      harmonic_success <- numeric(length(Qs))

      for (b in 1:B) {
        # Simulate a fresh trajectory for this replication.
        X <- numeric(n + burn)
        previous <- 0
        for (t in time) {
          if (case == 1) {
            alpha_t <- 0.7 * abs(sin(2 * pi * (t %% 7) / 7))
            innovation_mean <- 1
          } else {
            alpha_t <- 0.7
            # Exact positive-sine indicators; zeros at multiples of pi stay zero.
            r3 <- t %% 3
            r8 <- t %% 8
            innovation_mean <- 1 + 3 * (r3 > 0 & 2 * r3 < 3) +
                                    2 * (r8 > 0 & 2 * r8 < 8)
          }
          survivors <- rbinom(1, previous, alpha_t)
          if (case == 2) {
            p_t <- c(0.4, 0.5, 0.6, 0.7, 0.8)[((t - 1) %% 5) + 1]
            innovation <- rgeom(1, p_t)
          } else {
            innovation <- rpois(1, innovation_mean)
          }
          X[t] <- survivors + innovation
          previous <- X[t]
        }
        X <- tail(X, n)

        # Raw DFT periodogram: demean, omit zero and negative frequencies.
        power <- Mod(fft(X - mean(X)))^2 / n
        power <- power[2:(floor(n / 2) + 1)]

        # Automatic Fourier score for every candidate period.
        dft_strength <- numeric(max_Q)
        for (q in 2:max_Q) dft_strength[q] <- mean(power[harmonic_bins[[q]]])

        # Draw folds once; reuse them for every Q in this replication.
        foldid <- sample(rep(1:10, length.out = n))

        for (i in 1:length(Qs)) {
          Q <- Qs[i]
          keep <- column_period <= Q
          strength <- numeric(Q)
          if (var(X) > 0) {
            cv_fit <- glmnet::cv.glmnet(A[, keep, drop = FALSE], X,
              family = "gaussian", alpha = 1, intercept = TRUE,
              standardize = TRUE, type.measure = "mse", foldid = foldid,
              nlambda = 100, lambda.min.ratio = 1e-4, maxit = 1000000)
            if (cv_fit$glmnet.fit$jerr != 0)
              stop("glmnet convergence failure; do not discard this replication.")
            s <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
            s[abs(s) <= 1e-8] <- 0  # Numerical zero, not a significance threshold.
            for (q in 2:Q) strength[q] <- sum(abs(s[column_period[keep] == q]))
          }
          rpt_periods <- order(strength[2:Q], decreasing = TRUE) + 1
          rpt_periods <- head(rpt_periods[strength[rpt_periods] > 0], K)
          if (length(rpt_periods) == K && all(true_periods %in% rpt_periods))
            rpt_success[i] <- rpt_success[i] + 1

          dft_periods <- order(dft_strength[2:Q], decreasing = TRUE) + 1
          dft_periods <- head(dft_periods[dft_strength[dft_periods] > 0], K)
          if (length(dft_periods) == K && all(true_periods %in% dft_periods))
            harmonic_success[i] <- harmonic_success[i] + 1
        }
        if (b %% 50 == 0 || b == B)
          message("Case ", case, ", N=", n, ": ", b, "/", B)
      }

      for (i in 1:length(Qs)) {
        summary <- rbind(summary, data.frame(
          Case = case, N = n, Q = Qs[i], Replications = B,
          RPT_percent = 100 * rpt_success[i] / B,
          DFT_harmonic_percent = 100 * harmonic_success[i] / B,
          row.names = NULL))
      }
    }
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(outdir, "recovery_percentages.csv"), row.names = FALSE)
  print(summary, row.names = FALSE)
  message("Saved: ", file.path(outdir, "recovery_percentages.csv"))
  invisible(list(summary = summary))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  run_rpt_study(B = if (length(args)) as.integer(args[1]) else 500)
}
