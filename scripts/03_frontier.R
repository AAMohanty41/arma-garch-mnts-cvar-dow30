library(tidyverse)
library(quadprog)
library(yaml)
library(readr)
library(scales)
library(ggplot2)

# --- load inputs ---
cfg <- yaml::read_yaml("config/config.yml")
cap <- cfg$constraints$per_name_cap

obj <- readRDS("data/processed/mu_sigma.rds")
mu     <- as.numeric(obj$mu)         # monthly mean vector
Sigma0 <- obj$Sigma                  # monthly covariance
tickers <- obj$tickers
N <- length(mu)

# make covariance safely PD
eps <- max(1e-8, 1e-6 * mean(diag(Sigma0)))
Sigma <- Sigma0 + diag(eps, N)

# helpers -------------------------------------------------------------
vol_of <- function(w) sqrt(drop(t(w) %*% Sigma %*% w))
ret_of <- function(w) drop(sum(mu * w))

# Build constraint matrices for solve.QP: t(Amat) %*% w >= bvec
ones <- rep(1, N)
I    <- diag(N)

Amat_base <- cbind(ones,  I,  -I)                 # sum(w)=1; w>=0; -w>=-cap  (w<=cap)
bvec_base <- c(1,       rep(0, N), rep(-cap, N))
meq_base  <- 1                                    # first constraint is equality

# global min-variance (GMV), with bounds but no return target
gmv <- solve.QP(Dmat = 2*Sigma, dvec = rep(0, N),
                Amat = Amat_base, bvec = bvec_base, meq = meq_base)$solution

# extreme returns under the same bounds (approx LP with tiny quad term)
Dtiny <- diag(1e-8, N)
w_maxret <- solve.QP(Dmat = Dtiny, dvec =  mu,
                     Amat = Amat_base, bvec = bvec_base, meq = meq_base)$solution
w_minret <- solve.QP(Dmat = Dtiny, dvec = -mu,
                     Amat = Amat_base, bvec = bvec_base, meq = meq_base)$solution
r_max <- ret_of(w_maxret)
r_min <- ret_of(w_minret)

# grid of target returns across the feasible range
r_grid <- seq(r_min, r_max, length.out = 60)

# efficient frontier: minimize variance s.t. return >= r_t, bounds, and fully invested
solve_frontier_point <- function(r_t) {
  Amat <- cbind(ones, mu, I, -I)
  bvec <- c(1,     r_t, rep(0, N), rep(-cap, N))
  sol <- try(solve.QP(Dmat = 2*Sigma, dvec = rep(0, N),
                      Amat = Amat, bvec = bvec, meq = 1), silent = TRUE)
  if (inherits(sol, "try-error")) return(NULL)
  w <- sol$solution
  tibble(r_month = ret_of(w), s_month = vol_of(w)) |> mutate(type = "frontier")
}

frontier <- purrr::map(r_grid, solve_frontier_point) |> list_rbind() |> distinct()

# reference portfolios
w_ew <- rep(1/N, N)

refs <- tibble(
  label   = c("GMV", "Equal-Weight"),
  r_month = c(ret_of(gmv), ret_of(w_ew)),
  s_month = c(vol_of(gmv), vol_of(w_ew))
) |> mutate(type = "ref")

# combine and annualize for the table
out_tbl <- bind_rows(
  frontier,
  refs
) |>
  mutate(
    r_annual = r_month * 12,
    s_annual = s_month * sqrt(12)
  ) |>
  arrange(s_month, r_month)

# save table
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(out_tbl, "outputs/tables/efficient_frontier_points.csv")

# plot
p <- ggplot() +
  geom_path(data = frontier, aes(x = s_month, y = r_month), linewidth = 1) +
  geom_point(data = refs, aes(x = s_month, y = r_month), size = 2) +
  ggrepel::geom_text_repel(
    data = refs,
    aes(x = s_month, y = r_month, label = label),
    max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2, seed = 1
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    x = "Monthly volatility",
    y = "Monthly expected return",
    title = "Efficient Frontier (historical μ, Σ) • Long-only, w_i ≤ 10%"
  ) +
  theme_minimal(base_size = 12)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/efficient_frontier_past.png", p, width = 8, height = 5, dpi = 200)

cat("Frontier saved:\n",
    " - outputs/tables/efficient_frontier_points.csv\n",
    " - outputs/figures/efficient_frontier_past.png\n")

