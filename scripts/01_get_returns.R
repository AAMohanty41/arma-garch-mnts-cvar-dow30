library(tidyverse)
library(tidyquant)
library(lubridate)
library(readr)
library(yaml)

# 1) Load config
cfg <- yaml::read_yaml("config/config.yml")
tickers <- cfg$universe
start   <- as.Date(cfg$start_date)
end     <- as.Date(cfg$end_date)

# 2) Download adjusted prices
message("Downloading prices… chill.")
prices <- tq_get(tickers, from = start, to = end, get = "stock.prices", complete_cases = FALSE) |>
  select(symbol, date, adjusted)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
write_csv(prices, "data/raw/prices_adj.csv", na = "")

# 3) Monthly returns (simple + log)
monthly <- prices |>
  group_by(symbol) |>
  tq_transmute(
    select      = adjusted,
    mutate_fun  = periodReturn,
    period      = "monthly",
    type        = "arithmetic",
    col_rename  = "ret_simple"
  ) |>
  mutate(ret_log = log1p(ret_simple)) |>
  ungroup()

# 4) Align to complete panel (drop months with any NA across names)
wide <- monthly |>
  select(date, symbol, ret_log) |>
  pivot_wider(names_from = symbol, values_from = ret_log) |>
  arrange(date) |>
  drop_na()

# 5) Save long + wide
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(wide, "data/processed/returns_monthly_wide.csv", na = "")
long <- wide |>
  pivot_longer(-date, names_to = "symbol", values_to = "ret_log")
write_csv(long, "data/processed/returns_monthly_long.csv", na = "")






# === Daily log returns from adjusted prices ===
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

daily <- prices %>%
  arrange(symbol, date) %>%
  group_by(symbol) %>%
  mutate(
    ret_simple = adjusted/lag(adjusted) - 1,
    ret_log    = log1p(ret_simple)
  ) %>%
  ungroup()

# Align to a complete panel and drop dates with any missing ticker
wide_d <- daily %>%
  select(date, symbol, ret_log) %>%
  tidyr::pivot_wider(names_from = symbol, values_from = ret_log) %>%
  arrange(date) %>%
  tidyr::drop_na()

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(wide_d, "data/processed/returns_daily.rds")
readr::write_csv(wide_d, "data/processed/returns_daily.csv")
message("Saved daily returns -> data/processed/returns_daily.rds and .csv")




# 6) Tiny sanity report
kept    <- setdiff(names(wide), "date")
dropped <- setdiff(tickers, kept)
cat("\n=== RETURNS SUMMARY ===\n")
cat("Date range: ", as.character(min(wide$date)), " → ", as.character(max(wide$date)), "\n", sep = "")
cat("Months: ", nrow(wide), "\n", sep = "")
cat("Tickers kept (", length(kept), "): ", paste(kept, collapse = ", "), "\n", sep = "")
if (length(dropped)) cat("Dropped (insufficient history in window): ", paste(dropped, collapse = ", "), "\n", sep = "")
cat("Files:\n  - data/raw/prices_adj.csv\n  - data/processed/returns_monthly_wide.csv\n  - data/processed/returns_monthly_long.csv\n")
