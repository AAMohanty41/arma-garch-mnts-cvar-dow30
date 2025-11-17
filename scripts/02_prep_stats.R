library(tidyverse)
library(readr)

# load clean wide monthly log returns
wide <- read_csv("data/processed/returns_monthly_wide.csv", show_col_types = FALSE)

tickers <- setdiff(names(wide), "date")
X <- as.matrix(wide[tickers])        # T x N matrix of log returns
mu <- colMeans(X)                    # monthly mean vector
Sigma <- cov(X)                      # monthly covariance

# quick summary table
summary_tbl <- tibble(
  symbol = tickers,
  mean_monthly = mu,
  sd_monthly   = sqrt(diag(Sigma))
) |> arrange(desc(mean_monthly))

# save artifacts
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, "outputs/tables/summary_mu_sd.csv")
saveRDS(list(mu = mu, Sigma = Sigma, tickers = tickers), "data/processed/mu_sigma.rds")

cat("Stats saved:\n",
    " - outputs/tables/summary_mu_sd.csv\n",
    " - data/processed/mu_sigma.rds\n",
    "N =", length(tickers), "assets; months =", nrow(X), "\n")
