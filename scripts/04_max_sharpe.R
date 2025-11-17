library(tidyverse)
library(quadprog)
library(yaml)
library(readr)
library(ggplot2)
library(ggrepel)
library(scales)

# --- config & inputs ---
cfg <- yaml::read_yaml("config/config.yml")
cap <- cfg$constraints$per_name_cap
rf_ann <- if (!is.null(cfg$risk_free_annual)) cfg$risk_free_annual else 0
rf_month <- (1 + rf_ann)^(1/12) - 1

obj <- readRDS("data/processed/mu_sigma.rds")
mu      <- as.numeric(obj$mu)      # monthly mean vector
Sigma0  <- obj$Sigma               # monthly covariance
tickers <- obj$tickers
N <- length(mu)

# make covariance safely PD
eps <- max(1e-8, 1e-6 * mean(diag(Sigma0)))
Sigma <- Sigma0 + diag(eps, N)

clip_weights <- function(w, cap, tol = 1e-10) {
  w[abs(w) < tol] <- 0
  w <- pmax(0, pmin(w, cap))
  w / sum(w)
}

# apply to reference weights before saving
gmv  <- clip_weights(gmv,  cap)
w_ms <- clip_weights(w_ms, cap)
w_ew <- clip_weights(w_ew, cap)


# helpers
ret_of <- function(w) sum(mu * w)
vol_of <- function(w) sqrt(drop(t(w) %*% Sigma %*% w))

# constraints for solve.QP: t(Amat) %*% w >= bvec
ones <- rep(1, N); I <- diag(N)
Amat_base <- cbind(ones,  I,  -I)          # sum(w)=1; w>=0; -w>=-cap (w<=cap)
bvec_base <- c(1,       rep(0, N), rep(-cap, N))
meq_base  <- 1

# GMV, Max/Min return under bounds
gmv <- solve.QP(2*Sigma, rep(0, N), Amat_base, bvec_base, meq = meq_base)$solution
Dtiny <- diag(1e-8, N)
w_maxret <- solve.QP(Dtiny,  mu, Amat_base, bvec_base, meq = meq_base)$solution
w_minret <- solve.QP(Dtiny, -mu, Amat_base, bvec_base, meq = meq_base)$solution
r_max <- ret_of(w_maxret); r_min <- ret_of(w_minret)

# scan target-returns, keep weights
solve_point <- function(r_t) {
  Amat <- cbind(ones, mu, I, -I)
  bvec <- c(1,     r_t, rep(0, N), rep(-cap, N))
  sol <- try(solve.QP(2*Sigma, rep(0,N), Amat, bvec, meq = 1), silent = TRUE)
  if (inherits(sol, "try-error")) return(NULL)
  w <- pmax(0, pmin(sol$solution, cap))   # clip tiny numerical fuzz
  tibble(r = ret_of(w), s = vol_of(w),
         sharpe = ifelse(vol_of(w) > 0, (ret_of(w) - rf_month)/vol_of(w), NA_real_),
         w = list(setNames(as.numeric(w), tickers)))
}

r_grid <- seq(r_min, r_max, length.out = 80)
frontier_df <- purrr::map(r_grid, solve_point) |> bind_rows() |> distinct(r, s, .keep_all = TRUE)

# pick Max Sharpe
ms_idx <- which.max(frontier_df$sharpe)
w_ms <- frontier_df$w[[ms_idx]]
r_ms <- frontier_df$r[ms_idx]; s_ms <- frontier_df$s[ms_idx]; sh_ms <- frontier_df$sharpe[ms_idx]

# refs
w_ew <- rep(1/N, N)

# save weights
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(tibble(symbol = tickers, weight = w_ew), "outputs/tables/weights_equal_weight.csv")
write_csv(tibble(symbol = tickers, weight = gmv),  "outputs/tables/weights_gmv.csv")
write_csv(tibble(symbol = tickers, weight = w_ms), "outputs/tables/weights_max_sharpe.csv")

# save points (frontier + refs)
frontier_pts <- tibble(r_month = frontier_df$r, s_month = frontier_df$s, type = "frontier")
refs <- tibble(
  label   = c("GMV", "Equal-Weight", "Max Sharpe"),
  r_month = c(ret_of(gmv), ret_of(w_ew), r_ms),
  s_month = c(vol_of(gmv), vol_of(w_ew), s_ms),
  type = "ref"
)
write_csv(bind_rows(frontier_pts, refs), "outputs/tables/efficient_frontier_with_ms.csv")

# plot with Max Sharpe highlighted
# plot with legend explaining Max Sharpe
p <- ggplot() +
  geom_path(data = frontier_pts, aes(x = s_month, y = r_month, color = "Frontier"), linewidth = 1) +
  geom_point(data = refs, aes(x = s_month, y = r_month, color = label), size = 3) +
  scale_color_manual(
    name = "Legend",
    values = c("Frontier" = "black", "GMV" = "gray40", "Equal-Weight" = "steelblue", "Max Sharpe" = "red3"),
    breaks = c("Frontier", "GMV", "Equal-Weight", "Max Sharpe"),
    labels = c(
      "Frontier" = "Efficient Frontier",
      "GMV" = "Global Minimum Variance",
      "Equal-Weight" = "Equal-Weight Portfolio",
      "Max Sharpe" = "Max Sharpe (Highest Risk/Return Ratio)"
    )
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    x = "Monthly Volatility",
    y = "Monthly Expected Return",
    title = paste0("Efficient Frontier w/ Max Sharpe • rf = ",
                   percent(rf_month, accuracy = 0.01),
                   " per month; Long-only, w_i ≤ ", percent(cap))
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/efficient_frontier_past_with_MS_legend.png", p, width = 8, height = 5, dpi = 200)

cat("\nMax Sharpe (monthly rf =", sprintf("%.4f", rf_month), "):",
    "\n  r =", sprintf("%.4f", r_ms),
    " s =", sprintf("%.4f", s_ms),
    " SR =", sprintf("%.3f", sh_ms),
    "\nFiles:\n - outputs/tables/weights_max_sharpe.csv",
    "\n - outputs/tables/weights_gmv.csv",
    "\n - outputs/tables/weights_equal_weight.csv",
    "\n - outputs/tables/efficient_frontier_with_ms.csv",
    "\n - outputs/figures/efficient_frontier_past_with_MS.png\n")
