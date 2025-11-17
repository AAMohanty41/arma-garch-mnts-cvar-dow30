#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(tidyquant)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

message("=== 10_plot_jan2025_performance.R : Jan 21-day performance visuals ===")

# ------------------------------------------------------------
# 1) Load Jan 2025 summary table (from 09_compare_jan2025.R)
# ------------------------------------------------------------

summary_path <- "outputs/tables/jan2025_backtest_summary.csv"
if (!file.exists(summary_path)) {
  stop("Cannot find ", summary_path,
       " — run scripts/09_compare_jan2025.R first.")
}

jan_summary <- read.csv(summary_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 2) Barplot: predicted vs realized 21-day return
# ------------------------------------------------------------

df_long <- jan_summary %>%
  select(portfolio, pred_mean_21d, realized_21d) %>%
  pivot_longer(
    cols = c(pred_mean_21d, realized_21d),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    type = recode(
      type,
      pred_mean_21d = "Predicted mean (21d)",
      realized_21d  = "Realized return (21d)"
    )
  )

p_bar <- ggplot(df_long,
                aes(x = portfolio, y = value, fill = type)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    x = NULL,
    y = "21-day return",
    fill = NULL,
    title = "Predicted vs realized 21-day return\nHistorical vs CVaR Max-Sharpe portfolios"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("outputs/figures/jan2025_pred_vs_realized_returns.png",
       plot = p_bar,
       width = 7,
       height = 4,
       dpi   = 300)

message("Saved bar plot -> outputs/figures/jan2025_pred_vs_realized_returns.png")

# ------------------------------------------------------------
# 3) Time series: cumulative returns over the same 21-day window
# ------------------------------------------------------------

cfg <- yaml::read_yaml("config/config.yml")
tickers <- cfg$universe
last_train_date <- as.Date(cfg$end_date)  # e.g. "2024-12-31"

# ---- Load weights for historical Max-Sharpe ----

w_ms_path <- "outputs/tables/weights_max_sharpe.csv"
if (!file.exists(w_ms_path)) {
  stop("Cannot find ", w_ms_path, " — run 04_max_sharpe.R first.")
}

w_ms_df <- read.csv(w_ms_path, stringsAsFactors = FALSE)

# detect name/weight columns flexibly
name_col <- NULL
if ("ticker" %in% names(w_ms_df)) {
  name_col <- "ticker"
} else if ("symbol" %in% names(w_ms_df)) {
  name_col <- "symbol"
} else {
  stop("weights_max_sharpe.csv must have a 'ticker' or 'symbol' column; found: ",
       paste(names(w_ms_df), collapse = ", "))
}

weight_col <- NULL
if ("weight" %in% names(w_ms_df)) {
  weight_col <- "weight"
} else if ("w" %in% names(w_ms_df)) {
  weight_col <- "w"
} else {
  stop("weights_max_sharpe.csv must have a 'weight' or 'w' column; found: ",
       paste(names(w_ms_df), collapse = ", "))
}

w_ms_vec <- setNames(w_ms_df[[weight_col]], w_ms_df[[name_col]])

# ---- Load CVaR-frontier Max-Sharpe weights ----

frontier_path <- "outputs/tables/frontier_cvar_21d.csv"
if (!file.exists(frontier_path)) {
  stop("Cannot find ", frontier_path, " — run 07_frontier_cvar.R first.")
}

frontier <- read.csv(frontier_path, stringsAsFactors = FALSE)

if (!("Sharpe_sd" %in% names(frontier))) {
  stop("frontier_cvar_21d.csv has no 'Sharpe_sd' column.")
}

best_idx <- which.max(frontier$Sharpe_sd)

w_cols <- grep("^w_", names(frontier), value = TRUE)
if (length(w_cols) == 0) {
  stop("No w_* columns found in frontier_cvar_21d.csv.")
}

w_cvar_vec <- as.numeric(frontier[best_idx, w_cols])
names(w_cvar_vec) <- sub("^w_", "", w_cols)

# Enforce same ticker universe / order
common_tickers <- intersect(tickers, names(w_ms_vec))
common_tickers <- intersect(common_tickers, names(w_cvar_vec))
common_tickers <- sort(common_tickers)

w_ms_vec   <- w_ms_vec[common_tickers]
w_cvar_vec <- w_cvar_vec[common_tickers]

# ------------------------------------------------------------
# 4) Fetch the same 21 future trading days as in 09_compare_jan2025
# ------------------------------------------------------------

message("Downloading early-2025 prices for time series…")

prices_early <- tq_get(
  common_tickers,
  from = last_train_date - 15,
  to   = last_train_date + 60,
  get  = "stock.prices",
  complete_cases = FALSE
) %>%
  select(symbol, date, adjusted)

if (nrow(prices_early) == 0) {
  stop("No price data returned for early 2025. Check internet / ticker availability.")
}

daily_all <- prices_early %>%
  arrange(symbol, date) %>%
  group_by(symbol) %>%
  mutate(
    ret_log = log(adjusted / lag(adjusted))
  ) %>%
  ungroup()

future_dates <- sort(unique(daily_all$date[daily_all$date > last_train_date]))
h <- 21L
if (length(future_dates) < h) {
  stop("Fewer than ", h, " trading days after ", last_train_date,
       " in fetched data. Got: ", length(future_dates))
}

window_dates <- future_dates[1:h]

rets_win <- daily_all %>%
  filter(date %in% window_dates) %>%
  select(date, symbol, ret_log) %>%
  pivot_wider(names_from = symbol, values_from = ret_log) %>%
  arrange(date) %>%
  drop_na()

R_daily <- as.matrix(rets_win[, common_tickers])

# daily portfolio log returns
ret_ms_daily   <- as.numeric(R_daily %*% w_ms_vec)
ret_cvar_daily <- as.numeric(R_daily %*% w_cvar_vec)

# cumulative simple returns: exp(cumsum(log)) - 1
cum_ms   <- exp(cumsum(ret_ms_daily))   - 1
cum_cvar <- exp(cumsum(ret_cvar_daily)) - 1

ts_df <- tibble(
  date = rets_win$date,
  MS_historical    = cum_ms,
  MS_CVaR_frontier = cum_cvar
) %>%
  pivot_longer(-date, names_to = "portfolio", values_to = "cum_return")

p_ts <- ggplot(ts_df, aes(x = date, y = cum_return, color = portfolio)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    x = NULL,
    y = "Cumulative return over 21 days",
    color = NULL,
    title = "Cumulative performance over first 21 trading days after 2024-12-31"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("outputs/figures/jan2025_cum_returns_21d.png",
       plot = p_ts,
       width = 7,
       height = 4,
       dpi   = 300)

message("Saved time series plot -> outputs/figures/jan2025_cum_returns_21d.png")
message("=== Done: 10_plot_jan2025_performance.R ===")
