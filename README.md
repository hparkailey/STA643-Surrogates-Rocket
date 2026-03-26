# Turbine Blade Optimization — File Reference
**STA 643 Project 2**

---
 
## Pipeline Overview
 
```
generate_design.R
      |
      v
design_points.csv  -->  MATLAB simulator (simulator.p)
                                |
                                v
                   simulation_results.csv
                                |
                                v
                           fit_gp.R
                                |
                                v
                          gp_fits.RData  -->  [next: sequential BO]
```
 
---


## File Descriptions
 
### `generate_design.R`
**Purpose:** Generates the initial 100-point space-filling design using a
Maximin Latin Hypercube Design (MmLHD).
 
**Method:** Uses the `mymaximin_mmlhd()` function from HW3, which initializes
a Latin Hypercube Design and then runs a coordinate-swap optimization for
`T = 50000` iterations to maximize the minimum pairwise distance between
points. This gives both per-variable marginal coverage (the Latin Hypercube
property) and good joint space-filling (the maximin property).
 
A model-based design (e.g. MaxEnt, IMSPE) was deliberately avoided at this
stage because those require pre-specifying GP length-scales, which are unknown
prior to any data collection. MmLHD is more robust to this misspecification.
 
**Inputs:** None — design generated programmatically.
 
**Outputs:** `design_points.csv`


### `design_points.csv`
**Purpose:** The 100-point initial design in original (unscaled) physical
units, ready to be fed directly into the MATLAB simulator.
 
**Rows:** 100 (one per simulation run)


### `simulation_results.csv`
**Purpose:** The 100-point initial design with simulator outputs appended.
This is the primary input to `fit_gp.R`.

### `fit_gp.R`
**Purpose:** Fits two independent anisotropic Gaussian process (GP) surrogate
models to the simulation results and produces diagnostic plots.


 **Inputs:** `simulation_results.csv`
 
**Outputs:**
- `gp_fits.RData` — contains `gp_stress`, `gp_disp`, normalized design matrix
  `X`, log-transformed responses `y_stress` and `y_disp`, helper functions
  `normalize()`, `unnormalize()`, `predict_gp()`, and input bounds
- Diagnostic plots rendered to the active graphics device


## Budget Tracking
 
| Stage | Runs | Status |
|-------|------|--------|
| Initial MmLHD design | 100 | ✅ Complete |
| Sequential BO (constrained EI) | ~170 | ⬜ Next |
| Validation at optimum | ~30 | ⬜ Pending |
| **Total** | **300** | |
 
