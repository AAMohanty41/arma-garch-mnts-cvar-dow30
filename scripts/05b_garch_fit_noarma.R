# scripts/05_garch_fit.R  — ARMA-GARCH per asset -> eps + sigma
suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml", repos = "https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("forecast", quietly = TRUE)) install.packages("forecast", repos = "https://cloud.r-project.org")
  if (!requireNamespace("rugarch", quietly = TRUE)) install.packages("rugarch", repos = "https://cloud.r-project.org")
  library(yaml); library(data.table); library(forecast); library(rugarch)
})

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

load_cfg <- function() {
  cfg <- yaml::read_yaml("config/config.yml")
  if (is.null(cfg$universe)) stop("config.yml missing `universe:` list.")
  if (is.null(cfg$seed)) cfg$seed <- 42
  cfg
}

# Try to load saved returns; if missing, try to build them from 01 script or common objects
load_or_build_returns <- function(tickers) {
  ensure_dir("data/processed")
  candidates <- c("data/processed/returns_daily.rds",
                  "data/processed/returns_daily.csv",
                  "data/processed/returns.rds")
  found <- candidates[file.exists(candidates)]
  if (length(found)) {
    if (grepl("\\.rds$", found[1])) {
      DT <- as.data.table(readRDS(found[1]))
    } else {
      DT <- data.table::fread(found[1])
    }
  } else {
    # Try to source 01 to populate env
    if (file.exists("scripts/01_get_returns.R")) {
      try(source("scripts/01_get_returns.R"), silent = TRUE)
    }
    # Heuristic: look for an object that contains all tickers as columns
    objs <- mget(ls(envir = .GlobalEnv), envir = .GlobalEnv)
    DT <- NULL
    for (nm in names(objs)) {
      obj <- objs[[nm]]
      if (is.data.frame(obj) || is.matrix(obj)) {
        df <- as.data.frame(obj)
        cols <- intersect(tickers, colnames(df))
        if (length(cols) == length(tickers)) {
          DT <- as.data.table(df)
          break
        }
      }
    }
    if (is.null(DT)) {
      stop("Couldn't find daily returns. Run scripts/01_get_returns.R and make sure it creates a data.frame/matrix ",
           "with columns named by tickers. Then save it as data/processed/returns_daily.rds")
    }
    # Normalize Date if present and persist for next steps
    date_col <- intersect(names(DT), c("Date", "date"))
    if (length(date_col)) setnames(DT, date_col[1], "Date")
    if ("Date" %in% names(DT)) setorder(DT, Date)
    saveRDS(DT, "data/processed/returns_daily.rds")
  }
  
  date_col <- intersect(names(DT), c("Date", "date"))
  if (length(date_col)) setnames(DT, date_col[1], "Date")
  if ("Date" %in% names(DT)) setorder(DT, Date)
  
  cols <- intersect(tickers, names(DT))
  if (!length(cols)) stop("Loaded returns but none of the tickers are present as columns.")
  mask <- complete.cases(DT[, ..cols])
  
  list(Date = if ("Date" %in% names(DT)) DT$Date[mask] else NULL,
       X    = as.matrix(DT[mask, ..cols]),
       cols = cols)
}

fit_one <- function(x) {
  x <- as.numeric(x)
  a <- tryCatch(
    forecast::auto.arima(x, d = 0, max.p = 2, max.q = 2, seasonal = FALSE,
                         stationary = TRUE, ic = "aic", stepwise = FALSE, approximation = FALSE),
    error = function(e) NULL
  )
  p <- if (!is.null(a)) a$arma[1] else 0
  q <- if (!is.null(a)) a$arma[2] else 0
  
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  
  fit <- tryCatch(
    ugarchfit(spec = spec, data = x, solver = "hybrid", out.sample = 0, solver.control = list(trace = 0)),
    error = function(e) NULL
  )
  
  if (is.null(fit) || fit@fit$convergence != 0) {
    spec0 <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
      distribution.model = "norm"
    )
    fit <- ugarchfit(spec = spec0, data = x, solver = "hybrid", out.sample = 0, solver.control = list(trace = 0))
    p <- 0; q <- 0
  }
  
  z <- as.numeric(residuals(fit, standardize = TRUE))
  s <- as.numeric(sigma(fit))
  lb_p  <- tryCatch(Box.test(z,    lag = 12, type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  lb2_p <- tryCatch(Box.test(z^2,  lag = 12, type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  list(z = z, s = s, p = p, q = q, conv = fit@fit$convergence, lb = lb_p, lb2 = lb2_p)
}

# run ---------------------------------------------------------------------
cfg <- load_cfg()
set.seed(cfg$seed)
ensure_dir("data/processed"); ensure_dir("outputs/qa")

ret <- load_or_build_returns(cfg$universe)
X    <- ret$X; cols <- ret$cols; n <- nrow(X); k <- ncol(X)

eps   <- matrix(NA_real_, nrow = n, ncol = k, dimnames = list(NULL, cols))
sigma <- matrix(NA_real_, nrow = n, ncol = k, dimnames = list(NULL, cols))
qa    <- data.table(ticker = cols, p = NA_integer_, q = NA_integer_, convergence = NA_integer_,
                    lb_p = NA_real_, lb_p_sq = NA_real_)

for (j in seq_len(k)) {
  fitj <- fit_one(X[, j])
  eps[, j]   <- fitj$z
  sigma[, j] <- fitj$s
  qa$p[j]           <- fitj$p
  qa$q[j]           <- fitj$q
  qa$convergence[j] <- fitj$conv
  qa$lb_p[j]        <- fitj$lb
  qa$lb_p_sq[j]     <- fitj$lb2
}

if (!is.null(ret$Date)) {
  eps_dt   <- data.table(Date = ret$Date, eps)
  sigma_dt <- data.table(Date = ret$Date, sigma)
} else {
  eps_dt   <- as.data.table(eps)
  sigma_dt <- as.data.table(sigma)
}

saveRDS(eps_dt,   file = "data/processed/eps_daily_noarma.rds")
saveRDS(sigma_dt, file = "data/processed/sigma_daily_noarma.rds")
fwrite(qa,        file = "outputs/qa/garch_fit_summary_noarma.csv")

message("Saved: data/processed/eps_daily_noarma.rds")
message("Saved: data/processed/sigma_daily_noarma.rds")
message("Saved: outputs/qa/garch_fit_summary_noarma.csv")
