library(tidyverse)
library(readr)

## Inputs
# Monthly log returns (complete panel)
wide <- read_csv("data/processed/returns_monthly_wide.csv", show_col_types = FALSE)

tickers <- setdiff(names(wide), "date")
X <- as.matrix(wide[tickers])        # rows = months, cols = assets
mu <- colMeans(X)                    # monthly mean vector
Sigma <- cov(X)                      # monthly covariance

## Summary (mu, sigma)
summary_tbl <- tibble(
  symbol = tickers,
  mean_monthly = mu,
  sd_monthly   = sqrt(diag(Sigma))
) |> arrange(desc(mean_monthly))

## Outputs
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, "outputs/tables/summary_mu_sd.csv")
saveRDS(list(mu = mu, Sigma = Sigma, tickers = tickers), "data/processed/mu_sigma.rds")

cat("Stats saved:\n",
    " - outputs/tables/summary_mu_sd.csv\n",
    " - data/processed/mu_sigma.rds\n",
    "N =", length(tickers), "assets; months =", nrow(X), "\n")

