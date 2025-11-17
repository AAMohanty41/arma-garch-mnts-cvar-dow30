#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(CVXR)   # install.packages("CVXR") once if you don't have it
})

message("=== 07_frontier_cvar.R : CVaR frontier from MNTS sims ===")

# -------------------------------------------------------------------
# 1) Load config and scenarios
# -------------------------------------------------------------------

cfg <- yaml::read_yaml("config/config.yml")

alpha   <- if (!is.null(cfg$alpha))        cfg$alpha   else 0.95
w_cap   <- if (!is.null(cfg$w_cap))        cfg$w_cap   else 0.10
seed    <- if (!is.null(cfg$seed))         cfg$seed    else 42
long_only <- if (!is.null(cfg$long_only))  cfg$long_only else TRUE

set.seed(seed)

sim_path <- "data/processed/sims_mnts_21d.rds"
if (!file.exists(sim_path)) {
  stop("Cannot find ", sim_path, " — check data/processed/")
}

sims <- readRDS(sim_path)
R_mat <- as.matrix(sims)

S        <- nrow(R_mat)    # number of scenarios
n_assets <- ncol(R_mat)    # number of assets
tickers  <- colnames(R_mat)

if (is.null(tickers)) {
  stop("sims_mnts_21d.rds has no column names (tickers). Please fix upstream.")
}

message("Scenarios: ", S, "  |  Assets: ", n_assets)

# Sample mean per asset over 21-day horizon
mu_hat <- colMeans(R_mat)

# -------------------------------------------------------------------
# 2) Helper: compute stats for a given weight vector
# -------------------------------------------------------------------

portfolio_stats <- function(w, R_mat, alpha = 0.95) {
  # w: numeric vector of weights (length = ncol(R_mat))
  port_ret <- as.numeric(R_mat %*% w)  # 21-day returns per scenario
  
  mean_ret <- mean(port_ret)
  sd_ret   <- stats::sd(port_ret)
  
  # Losses = -return (positive = bad)
  losses <- -port_ret
  
  var_alpha  <- unname(stats::quantile(losses, alpha, type = 7))
  cvar_alpha <- mean(losses[losses >= var_alpha])
  
  list(
    mean = mean_ret,
    sd   = sd_ret,
    var  = var_alpha,
    cvar = cvar_alpha
  )
}

# -------------------------------------------------------------------
# 3) Optimization: minimize CVaR for a target mean
#    Rockafellar–Uryasev formulation with CVXR
# -------------------------------------------------------------------

solve_cvar_for_target <- function(mu_target,
                                  R_mat,
                                  mu_hat,
                                  alpha,
                                  w_cap,
                                  long_only = TRUE) {
  S        <- nrow(R_mat)
  n_assets <- ncol(R_mat)
  
  w     <- CVXR::Variable(n_assets)  # portfolio weights
  t_var <- CVXR::Variable(1)         # VaR threshold
  u     <- CVXR::Variable(S)         # excess losses
  
  # Portfolio losses for each scenario: ℓ_s = -R_s * w
  losses <- - (R_mat %*% w)
  
  objective <- t_var + (1 / ((1 - alpha) * S)) * sum(u)
  
  constr <- list(
    sum(w) == 1,
    u >= 0,
    u >= losses - t_var,
    t(mu_hat) %*% w >= mu_target
  )
  
  if (long_only) {
    constr <- c(constr, list(
      w >= 0,
      w <= w_cap
    ))
  } else {
    # If you ever allow shorting, you can relax the lower bound
    constr <- c(constr, list(
      w >= -w_cap,
      w <=  w_cap
    ))
  }
  
  prob   <- CVXR::Problem(CVXR::Minimize(objective), constr)
  result <- CVXR::solve(prob, solver = "ECOS")
  
  if (result$status != "optimal") {
    message("Target ", signif(mu_target, 4),
            " infeasible / not optimal (status: ", result$status, ")")
    return(NULL)
  }
  
  w_opt <- as.numeric(result$getValue(w))
  names(w_opt) <- colnames(R_mat)
  
  stats <- portfolio_stats(w_opt, R_mat, alpha)
  
  list(
    mu_target = mu_target,
    weights   = w_opt,
    stats     = stats
  )
}

# -------------------------------------------------------------------
# 4) Choose a sensible grid of target means
# -------------------------------------------------------------------

# Use equal-weight portfolio as a reference point
w_eq      <- rep(1 / n_assets, n_assets)
eq_stats  <- portfolio_stats(w_eq, R_mat, alpha)
mu_center <- eq_stats$mean

# Go from ~half to ~1.5x that mean; keep order increasing
mu_low  <- mu_center * 0.5
mu_high <- mu_center * 1.5
if (mu_low > mu_high) {
  tmp     <- mu_low
  mu_low  <- mu_high
  mu_high <- tmp
}

mu_targets <- seq(mu_low, mu_high, length.out = 25)

message("Target mean range (21d): [",
        signif(mu_low, 5), ", ",
        signif(mu_high, 5), "] with ",
        length(mu_targets), " grid points")

# -------------------------------------------------------------------
# 5) Solve CVaR problem across the grid
# -------------------------------------------------------------------

solutions <- list()
for (mu_t in mu_targets) {
  message("Solving CVaR frontier point for target mean = ",
          signif(mu_t, 5), " ...")
  sol <- solve_cvar_for_target(mu_t, R_mat, mu_hat, alpha, w_cap, long_only)
  if (!is.null(sol)) {
    solutions[[length(solutions) + 1]] <- sol
  }
}

if (length(solutions) == 0) {
  stop("No feasible CVaR frontier points found. Check constraints / target range.")
}

# -------------------------------------------------------------------
# 6) Build output table with metrics + weights
# -------------------------------------------------------------------

rows <- lapply(solutions, function(sol) {
  w     <- sol$weights
  stats <- sol$stats
  
  sharpe_sd   <- if (stats$sd   > 0) stats$mean / stats$sd   else NA_real_
  sharpe_cvar <- if (stats$cvar > 0) stats$mean / stats$cvar else NA_real_
  
  row <- data.frame(
    target_mu  = sol$mu_target,
    mean_21d   = stats$mean,
    sd_21d     = stats$sd,
    VaR95_loss = stats$var,
    CVaR95_loss = stats$cvar,
    Sharpe_sd   = sharpe_sd,
    Sharpe_CVaR = sharpe_cvar,
    stringsAsFactors = FALSE
  )
  
  # Attach weights as columns w_TICKER
  for (nm in names(w)) {
    row[[paste0("w_", nm)]] <- w[[nm]]
  }
  
  row
})

frontier_cvar <- do.call(rbind, rows)

# -------------------------------------------------------------------
# 7) Save to CSV and report max Sharpe
# -------------------------------------------------------------------

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

out_path <- "outputs/tables/frontier_cvar_21d.csv"
utils::write.csv(frontier_cvar, out_path, row.names = FALSE)

message("Written CVaR frontier table to: ", out_path)

# Find max Sharpe (sd-based) on this frontier
best_idx <- which.max(frontier_cvar$Sharpe_sd)
if (length(best_idx) == 1 && is.finite(frontier_cvar$Sharpe_sd[best_idx])) {
  message("Max Sharpe (sd-based) on CVaR frontier at row ", best_idx, ":")
  message("  mean_21d   = ", signif(frontier_cvar$mean_21d[best_idx], 6))
  message("  sd_21d     = ", signif(frontier_cvar$sd_21d[best_idx], 6))
  message("  CVaR95_loss= ", signif(frontier_cvar$CVaR95_loss[best_idx], 6))
} else {
  message("Could not identify a unique max Sharpe point (check data).")
}

message("=== Done: CVaR frontier from MNTS 21d sims ===")
