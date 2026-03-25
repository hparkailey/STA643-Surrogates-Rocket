# =============================================================================
# generate_design.R
# Generates MmLHD initial design points for the turbine blade simulator
#
# OUTPUT: design_points.csv — rows are runs, columns are the 6 simulator inputs
#         in their ORIGINAL (unscaled) units, ready to feed into MATLAB
# =============================================================================

library(plgp)

# =============================================================================
# 1.MmLHD FUNCTION (from HW3.Rmd)
# =============================================================================

mymaximin_mmlhd <- function(n, m, T = 100000) {
  
  # STEP 1: Initialize as LHD
  l <- (-(n - 1)/2):((n - 1)/2)
  L <- matrix(NA, nrow = n, ncol = m)
  for (j in 1:m) L[, j] <- sample(l, n)
  X <- (L + (n - 1)/2 + 0.5) / n  # points in [0, 1]^m
  
  d   <- distance(X)
  d   <- d[upper.tri(d)]
  md  <- min(d)
  
  # STEPS 2-4: Swap coordinates to maximize minimum distance
  for (t in 1:T) {
    rows <- sample(1:n, 2)
    col  <- sample(1:m, 1)
    
    # Propose swap
    val1 <- X[rows[1], col]
    val2 <- X[rows[2], col]
    X[rows[1], col] <- val2
    X[rows[2], col] <- val1
    
    d_mat   <- distance(X)
    d_mat   <- d_mat[upper.tri(d_mat)]
    mdprime <- min(d_mat)
    
    # Accept if improved, otherwise revert
    if (mdprime > md) {
      md <- mdprime
    } else {
      X[rows[1], col] <- val1
      X[rows[2], col] <- val2
    }
  }
  
  return(X)
}

# =============================================================================
# 2. DESIGN PARAMETERS
# =============================================================================

n_runs <- 100   # initial budget  (~1/3 of total 300)
d      <- 6     # number of input variables
T_iter <- 50000 # swap iterations (same as your HW3 Q3)

set.seed(1)

# =============================================================================
# 3. INPUT VARIABLE RANGES 
# =============================================================================
#   x1: Young's Modulus          [200e9,  300e9]  Pa
#   x2: Poisson's Ratio          [0.1,    0.49]   (unitless)
#   x3: Coeff. Thermal Expansion [5e-6,   15e-6]  K^-1
#   x4: Thermal Conductivity     [5,      15]     W/m/K
#   x5: Internal Cooling Temp    [50,     350]    deg C
#   x6: Pressure Load on Suction [1e5,    4.8e5]  Pa

lower <- c(200e9,  0.10,  5e-6,  5,   50,   1.0e5)
upper <- c(300e9,  0.49,  15e-6, 15,  350,  4.8e5)

var_names <- c("ymod", "prat", "cte", "therm", "ctemp", "press")

# =============================================================================
# 4. GENERATE MmLHD IN [0,1]^6 THEN SCALE TO ORIGINAL UNITS
# =============================================================================

cat("Generating", n_runs, "-point MmLHD in d =", d, "dimensions...\n")
cat("(T =", T_iter, "swap iterations)\n")

X_unit <- mymaximin_mmlhd(n = n_runs, m = d, T = T_iter)

# Scale each column from [0,1] to [lower_l, upper_l]
X_scaled <- X_unit
for (l in 1:d) {
  X_scaled[, l] <- lower[l] + X_unit[, l] * (upper[l] - lower[l])
}

colnames(X_scaled) <- var_names

# =============================================================================
# 5. QUICK DIAGNOSTICS
# =============================================================================

# Maximin criterion on unit design (useful for comparing vs random LHD)
get_md <- function(X) {
  dst <- distance(X)
  min(dst[upper.tri(dst)])
}

cat("\nMaximin criterion (unit scale):", round(get_md(X_unit), 4), "\n")
cat("Design dimensions:", nrow(X_scaled), "runs x", ncol(X_scaled), "variables\n")
cat("\nFirst 5 rows (original units):\n")
print(round(X_scaled[1:5, ], 6))

# Optional: verify marginal coverage looks uniform (each variable)
cat("\nColumn ranges (should match input bounds):\n")
for (l in 1:d) {
  cat(sprintf("  %s: [%.4g, %.4g]\n", var_names[l], 
              min(X_scaled[, l]), max(X_scaled[, l])))
}

# =============================================================================
# 6. SAVE TO CSV FOR MATLAB
# =============================================================================
# Each row = one simulator run
# Columns: ymod, prat, cte, therm, ctemp, press (in original units)

write.csv(X_scaled, "design_points.csv", row.names = FALSE)
cat("\nSaved to design_points.csv\n")

# =============================================================================
# 7. HOW TO READ THIS IN MATLAB
# =============================================================================
# In MATLAB, run:
#
#   design = readtable('design_points.csv');
#   results = zeros(height(design), 2);  % cols: stress, displacement
#
#   for i = 1:height(design)
#       ymod  = design.ymod(i);
#       prat  = design.prat(i);
#       cte   = design.cte(i);
#       therm = design.therm(i);
#       ctemp = design.ctemp(i);
#       press = design.press(i);
#       [stress, displ] = simulator(ymod, prat, cte, therm, ctemp, press);
#       results(i, :) = [stress, displ];
#       fprintf('Run %d/%d complete\n', i, height(design));
#   end
#
#   % Save outputs
#   results_table = array2table(results, 'VariableNames', {'stress','displacement'});
#   writetable([design, results_table], 'simulation_results.csv');
