# ARMA–GARCH + MNTS Mean–CVaR Frontier (Dow 30)

This project builds and compares two efficient frontiers for portfolios on the Dow 30:

- A historical Markowitz mean–variance frontier using 5 years of monthly returns.
- A model-based mean–CVaR frontier using ARMA–GARCH volatility modeling and a multivariate Normal Tempered Stable (MNTS) distribution with Monte Carlo simulation.

The goal is to study how a tail-risk-aware, CVaR-optimized portfolio compares to the historical Max-Sharpe portfolio, and to check how well the model predicts out-of-sample performance.

## Universe & horizon

- Universe: Dow 30 constituents
- Data: daily prices over a 5-year window ending 2024-12-31
- Horizon: 21 trading days
- Scenarios: 20,000 simulated 21-day paths

Configuration (dates, tickers, etc.) is stored in `config/config.yml`.

## Pipeline overview

The scripts are meant to be run in roughly numeric order:

1. `scripts/01_get_returns.R`  
   - Downloads adjusted prices for the Dow 30 using `tidyquant`.  
   - Computes daily and monthly log returns.  
   - Saves `data/raw/prices_adj.csv` and `data/processed/returns_*`.

2. `scripts/02_prep_stats.R`  
   - Computes sample means/standard deviations and other summary stats for monthly returns.  
   - Prepares inputs for Markowitz optimization.

3. `scripts/03_frontier.R` 
   - Builds the historical Markowitz mean–variance efficient frontier (long-only, 0–10% per-name cap).  
   - Exports `outputs/tables/efficient_frontier_points.csv`.

4. `scripts/04_max_sharpe.R`  
   - Identifies Equal-Weight, Global Minimum Variance (GMV), and Max-Sharpe portfolios on the historical frontier.  
   - Saves their weights in `outputs/tables/weights_*.csv` and a summary in `outputs/tables/summary_mu_sd.csv`.  
   - Produces the plot `outputs/figures/efficient_frontier_past_with_MS_legend.png`.

5. `scripts/05_garch_fit.R`  
   - Fits ARMA–GARCH models to each stock’s daily returns.  
   - Saves standardized residuals and conditional volatilities for use in the MNTS fit.

6. `scripts/06_mnts_fit_sim.R`  
   - Fits a multivariate Normal Tempered Stable (MNTS) distribution to the standardized residuals.  
   - Simulates 20,000 joint 21-day return paths for the Dow 30.  
   - Exports simulation objects under `data/processed/` and QA summaries under `outputs/qa/`.

7. `scripts/07_frontier_cvar.R`**  
   - Uses the simulated 21-day returns to construct a mean–CVaR efficient frontier for 25 target-return levels.  
   - Optimizes portfolios using a CVaR objective (via `CVXR`).  
   - Saves `outputs/tables/frontier_cvar_21d.csv` and reports the Max-Sharpe point on this frontier.

8. `scripts/08_plot_cvar_frontier.R`  
   - Plots 21-day expected return vs 95% CVaR (loss) for the model-based frontier.  
   - Highlights the Max-Sharpe portfolio.  
   - Exports `outputs/figures/frontier_cvar_21d.png`.

9. `scripts/09_compare_jan2025.R`  
   - Treats 2024-12-31 as the end of the training window.  
   - Fetches the first 21 trading days after this date (early 2025).  
   - Computes realized 21-day returns for:  
     - the historical Max-Sharpe portfolio, and  
     - the CVaR-optimized Max-Sharpe portfolio.  
   - Compares these realized returns to the simulated distribution (percentiles, VaR95 breaches).  
   - Writes a compact summary to `outputs/tables/jan2025_backtest_summary.csv`.

## Outputs

Key visual outputs:

- `outputs/figures/efficient_frontier_past_with_MS_legend.png`  
  Historical Markowitz efficient frontier with Equal-Weight, GMV, and Max-Sharpe portfolios.

- `outputs/figures/frontier_cvar_21d.png`  
  Model-based 21-day mean–CVaR frontier from ARMA–GARCH + MNTS simulations, with the Max-Sharpe portfolio highlighted.

Data, large simulation files, and QA tables are intentionally excluded from version control to keep the repository lightweight.
