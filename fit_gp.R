# =============================================================================
# fit_gp.R
# Fits two anisotropic GP surrogates to turbine blade simulation results:
#   GP1: max stress       (objective to minimize)
#   GP2: max displacement (constraint: must be < 1.3e-3)
#
# Both outputs are log-transformed before fitting (raw scale caused skewed
# residuals and optimizer boundary issues in earlier attempts).
#
# INPUTS:  simulation_results.csv  (from MATLAB simulator)
# OUTPUTS: gp_fits.RData           (fitted GP objects for use in BO)
# =============================================================================

library(plgp)   # distance()
library(hetGP)  # mleHomGP(), predict()

# =============================================================================
# 1. LOAD & INSPECT DATA
# =============================================================================

dat <- read.csv("simulation_results.csv")
cat("Loaded", nrow(dat), "simulation runs\n")
cat("Columns:", paste(names(dat), collapse = ", "), "\n")

# Separate inputs and outputs
input_names <- c("ymod", "prat", "cte", "therm", "ctemp", "press")
X_raw       <- as.matrix(dat[, input_names])

# BUG FIX 1: constraint check must use raw displacement, not log-transformed
cat("\nRaw stress      — min:", min(dat$stress),       "max:", max(dat$stress),       "\n")
cat("Raw displacement — min:", min(dat$displacement), "max:", max(dat$displacement), "\n")
cat("Runs violating displacement constraint (> 1.3e-3):",
    sum(dat$displacement > 1.3e-3), "/", nrow(dat), "\n")

# Log-transform outputs: fixes skewed residuals and optimizer boundary issues
# Both are strictly positive physical quantities so log is well-defined
y_stress <- log(dat$stress)
y_disp   <- log(dat$displacement)

cat("\nLog stress      — min:", round(min(y_stress), 3), "max:", round(max(y_stress), 3), "\n")
cat("Log displacement — min:", round(min(y_disp),   3), "max:", round(max(y_disp),   3), "\n")

# =============================================================================
# 2. NORMALIZE INPUTS TO [0,1]  (GP fitting must be on unit scale)
# =============================================================================

lower <- c(200e9, 0.10, 5e-6,  5,   50,  1.0e5)
upper <- c(300e9, 0.49, 15e-6, 15, 350,  4.8e5)

normalize <- function(X_raw, lower, upper) {
  X_unit <- X_raw
  for (l in 1:ncol(X_raw)) {
    X_unit[, l] <- (X_raw[, l] - lower[l]) / (upper[l] - lower[l])
  }
  return(X_unit)
}

unnormalize <- function(X_unit, lower, upper) {
  X_raw <- X_unit
  for (l in 1:ncol(X_unit)) {
    X_raw[, l] <- lower[l] + X_unit[, l] * (upper[l] - lower[l])
  }
  return(X_raw)
}

X <- normalize(X_raw, lower, upper)
cat("\nInputs normalized to [0,1]. Column ranges:\n")
for (l in 1:ncol(X)) cat(sprintf("  x%d: [%.3f, %.3f]\n", l, min(X[, l]), max(X[, l])))

# =============================================================================
# 3. GP FITTING & PREDICTION FUNCTIONS (hetGP)
#    Kernel: Gaussian (= squared exponential), anisotropic
#    Estimates: theta_1,...,theta_6 (per-dimension length-scales), nu_hat, g
#    NOTE: lower/upper here are bounds on theta hyperparameters,
#          NOT on the input variable ranges
# =============================================================================

fit_gp_mle <- function(X, y) {
  out <- mleHomGP(X       = X,
                  Z       = y,
                  lower   = rep(1e-3, ncol(X)),  # lower bounds on theta
                  upper   = rep(20,   ncol(X)),  # upper bounds on theta
                  covtype = "Gaussian")           # Gaussian = SE kernel
  return(out)
}

predict_gp <- function(gp, X_new) {
  pred <- predict(gp, x = X_new)
  return(list(
    mean     = pred$mean,
    variance = pred$sd2,
    sd       = sqrt(pred$sd2)
  ))
}

# =============================================================================
# 4. LOO-CV DIAGNOSTIC FUNCTION
# Note: hetGP does not expose Ki directly so we use in-sample predictions
# at training points as a proxy. The nugget g acts as regularization so
# in-sample predictions are not trivially perfect (unlike a noiseless GP).
# =============================================================================

loo_cv_hetgp <- function(gp, y) {
  # Predict at training locations (gp$X0 stores unique design points)
  pred      <- predict(gp, x = gp$X0)
  residuals <- y - pred$mean
  std_resid <- residuals / sqrt(pred$sd2)
  rmse      <- sqrt(mean(residuals^2))
  coverage  <- mean(abs(std_resid) < qnorm(0.95))  # should be ~0.90
  return(list(
    residuals   = residuals,
    std_resid   = std_resid,
    rmse        = rmse,
    coverage_90 = coverage,
    loo_mean    = pred$mean,
    loo_sd      = sqrt(pred$sd2)
  ))
}

# =============================================================================
# 5. FIT GP1 (LOG STRESS) AND GP2 (LOG DISPLACEMENT)
# =============================================================================

cat("\n--- Fitting GP1: Log Max Stress ---\n")
set.seed(643)
gp_stress <- fit_gp_mle(X, y_stress)
cat("Estimated length-scales (theta):\n")
for (l in 1:6) cat(sprintf("  theta_%d (%s): %.4f\n", l, input_names[l], gp_stress$theta[l]))
cat("mu (beta0) =",    round(gp_stress$beta0,  4),
    " tau2 (nu_hat) =", round(gp_stress$nu_hat, 6),
    " g =",             round(gp_stress$g, 8), "\n")

cat("\n--- Fitting GP2: Log Max Displacement ---\n")
set.seed(643)
gp_disp <- fit_gp_mle(X, y_disp)
cat("Estimated length-scales (theta):\n")
for (l in 1:6) cat(sprintf("  theta_%d (%s): %.4f\n", l, input_names[l], gp_disp$theta[l]))
cat("mu (beta0) =",    round(gp_disp$beta0,  6),
    " tau2 (nu_hat) =", round(gp_disp$nu_hat, 8),
    " g =",             round(gp_disp$g, 8), "\n")

# =============================================================================
# 6. MODEL DIAGNOSTICS
# =============================================================================

cat("\n--- Diagnostics ---\n")

# BUG FIX 4: pass the correct log-transformed vectors consistently
loo_s <- loo_cv_hetgp(gp_stress, y_stress)
loo_d <- loo_cv_hetgp(gp_disp,   y_disp)

cat(sprintf("GP1 (Log Stress):       RMSE = %.4f  |  90%% coverage = %.2f (target 0.90)\n",
            loo_s$rmse, loo_s$coverage_90))
cat(sprintf("GP2 (Log Displacement): RMSE = %.6f  |  90%% coverage = %.2f (target 0.90)\n",
            loo_d$rmse, loo_d$coverage_90))

# =============================================================================
# 7. DIAGNOSTIC PLOTS
# =============================================================================

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# --- GP1: Predicted vs Actual ---
# BUG FIX 2: gp_stress$y does not exist in hetGP objects — use y_stress directly
plot(y_stress, loo_s$loo_mean,
     xlab = "Actual Log Stress", ylab = "Predicted Log Stress",
     main = "GP1 (Stress): Predicted vs Actual", pch = 19, col = "#2b6ca8aa")
abline(0, 1, col = "red", lwd = 2)

# --- GP1: Standardized residuals (should look N(0,1)) ---
hist(loo_s$std_resid, breaks = 20, probability = TRUE,
     main = "GP1: Standardized Residuals",
     xlab = "Standardized Residual", col = "lightblue")
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

# --- GP1: Length-scale importance ---
barplot(gp_stress$theta,
        names.arg = input_names,
        main = "GP1 (Stress): Length-scales\n(larger = less important)",
        ylab = expression(theta[l]),
        col  = "#2b6ca8")

# --- GP2: Predicted vs Actual ---
# BUG FIX 2: gp_disp$y does not exist in hetGP objects — use y_disp directly
plot(y_disp, loo_d$loo_mean,
     xlab = "Actual Log Displacement", ylab = "Predicted Log Displacement",
     main = "GP2 (Displacement): Predicted vs Actual", pch = 19, col = "#a83232aa")
abline(0, 1, col = "red", lwd = 2)
# BUG FIX 3: constraint threshold must be on log scale since plot is log scale
abline(h = log(1.3e-3), col = "orange", lty = 2, lwd = 1.5)
legend("topleft", legend = "log(d*) = log(1.3e-3)", col = "orange", lty = 2, bty = "n")

# --- GP2: Standardized residuals ---
hist(loo_d$std_resid, breaks = 20, probability = TRUE,
     main = "GP2: Standardized Residuals",
     xlab = "Standardized Residual", col = "lightcoral")
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

# --- GP2: Length-scale importance ---
barplot(gp_disp$theta,
        names.arg = input_names,
        main = "GP2 (Displacement): Length-scales\n(larger = less important)",
        ylab = expression(theta[l]),
        col  = "#a83232")

# =============================================================================
# 8. VARIABLE IMPORTANCE SUMMARY
# =============================================================================

cat("\n--- Variable Importance (smaller theta = more important) ---\n")

cat("\nGP1 (Log Stress) ranking:\n")
rank_s <- order(gp_stress$theta)
for (i in 1:6) cat(sprintf("  %d. %s  (theta = %.4f)\n",
                            i, input_names[rank_s[i]], gp_stress$theta[rank_s[i]]))

cat("\nGP2 (Log Displacement) ranking:\n")
rank_d <- order(gp_disp$theta)
for (i in 1:6) cat(sprintf("  %d. %s  (theta = %.4f)\n",
                            i, input_names[rank_d[i]], gp_disp$theta[rank_d[i]]))

# =============================================================================
# 9. SAVE FITTED GP OBJECTS FOR USE IN SEQUENTIAL BO
# =============================================================================
# Note: y_stress and y_disp saved here are log-transformed.
# When back-transforming predictions in BO script, use exp() on predicted means.

save(gp_stress, gp_disp,
     X, y_stress, y_disp,
     lower, upper, input_names,
     normalize, unnormalize, predict_gp,
     file = "gp_fits.RData")

cat("\nSaved gp_stress, gp_disp to gp_fits.RData\n")
cat("Load in next script with: load('gp_fits.RData')\n")
cat("Note: y_stress and y_disp are on log scale — use exp() to back-transform.\n")
