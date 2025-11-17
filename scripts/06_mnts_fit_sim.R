# 06_mnts_fit_sim.R
# Fit MNTS with temStaR on daily returns and simulate 2-day horizon

suppressPackageStartupMessages({
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  if (!requireNamespace("temStaR", quietly = TRUE)) remotes::install_github("aaron9011/temStaR-v0.90")
  
  # CDF generic lives in spatstat.core
  for (p in c("spatstat.geom", "spatstat.random", "spatstat.explore", "spatstat.core", "data.table", "readr", "yaml")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  }
  
  library(temStaR)
  library(spatstat.core)
  library(data.table)
  library(readr)
  library(yaml)
})







`%||%` <- function(x, y) if (is.null(x)) y else x

# --- load config + ensure dirs ---
cfg  <- yaml::read_yaml("config/config.yml")
dirs <- c("data/processed", "outputs/qa", "outputs/tables", "outputs/sims")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)






# --- load data robustly ---
rd <- readRDS("data/processed/returns_daily.rds")

# coerce to plain data.frame and strip any date columns
rd <- as.data.frame(rd)
rd$Date <- NULL
rd$date <- NULL

# keep only numeric columns (drop anything weird)
num_cols <- sapply(rd, is.numeric)
rd <- rd[, num_cols, drop = FALSE]

# force to numeric matrix and drop any non-finite rows
X <- data.matrix(rd)
good <- stats::complete.cases(X) & apply(is.finite(X), 1, all)
X <- X[good, , drop = FALSE]

tickers <- colnames(X)
n <- nrow(X); k <- ncol(X)
stopifnot(n > 100, k >= 2)






# --- fit MNTS with temStaR ---
# This estimates mu, sigma, alpha, theta, beta, and Rho in one go.
message("Fitting MNTS on daily returns... this can take a minute.")
st <- temStaR::fitmnts(returndata = X, n = k)

# --- persist the fitted structure ---
saveRDS(st, file = "data/processed/mnts_fit.rds")

# --- QA summaries ---
# 1) global params
qa_global <- data.table(
  param = c("alpha", "theta"),
  value = c(st$alpha, st$theta)
)

# 2) per-asset marginals
qa_assets <- data.table(
  symbol = tickers,
  mu     = as.numeric(st$mu),
  sigma  = as.numeric(st$sigma),
  beta   = as.numeric(st$beta)
)

# 3) correlation matrix of stdNTS (epsilon)
qa_rho <- as.data.table(st$Rho)
setnames(qa_rho, names(qa_rho), tickers)
qa_rho <- cbind(symbol = tickers, qa_rho)

fwrite(qa_global, "outputs/qa/mnts_fit_summary_global.csv")
fwrite(qa_assets, "outputs/qa/mnts_fit_summary_assets.csv")
fwrite(qa_rho,    "outputs/qa/mnts_fit_rho.csv")

message("Saved: outputs/qa/mnts_fit_summary_global.csv")
message("Saved: outputs/qa/mnts_fit_summary_assets.csv")
message("Saved: outputs/qa/mnts_fit_rho.csv")

# --- simulate 2-day horizon ---
set.seed(cfg$seed %||% 42)
nsim <- cfg$scenarios %||% 20000
H    <- cfg$horizon_days %||% 2

# Draw H independent MNTS days and sum them
sim_list <- vector("list", H)
for (h in seq_len(H)) {
  sim_list[[h]] <- temStaR::rmnts(st, numofsample = nsim)
}
sim_agg <- Reduce(`+`, sim_list)    # nsim x k
colnames(sim_agg) <- tickers

# Save both RDS and CSV (CSV is large; keep it optional if needed)
saveRDS(sim_agg, "data/processed/sims_mnts_2d.rds")
fwrite(as.data.table(sim_agg), "outputs/tables/sims_mnts_2d.csv")

message("Saved: data/processed/mnts_fit.rds")
message("Saved: data/processed/sims_mnts_2d.rds")
message("Saved: outputs/tables/sims_mnts_2d.csv")

# --- tiny print to console so you see it's sane ---
cat(sprintf("MNTS fit: alpha=%.4f  theta=%.4f | k=%d assets | nsim=%d | horizon=%d days\n",
            st$alpha, st$theta, k, nsim, H))
