# AGENTS.md

## Proactive planning and follow-through

- When we work through a longer plan or todo list, track what is
  completed, what remains, and any blockers.
- Whenever you finish a meaningful part of the work, recommend the
  next concrete task and briefly explain why it should come next.
- If several next steps are plausible, suggest up to three in
  priority order and identify your recommendation.
- Make suggestions specific to our current project. Avoid generic
  offers such as "Let me know if you need anything else."
- Continue through already authorized steps without repeatedly
  asking whether to proceed.
- Present ideas that expand the agreed scope as suggestions.
- When the agreed plan is complete, say so. Do not invent extra work.

---

## Project Overview

This repository replicates the **NIST StRD (Statistical Reference Datasets)** numerical accuracy benchmark in Python, modeled after the [Stata 17 NIST StRD certification report](https://www.stata.com/support/cert/nist17/).

- **NIST StRD Data Archive**: `https://www.itl.nist.gov/div898/strd/general/dataarchive.html`
- **Reference Benchmark Report**: `https://www.stata.com/support/cert/nist17/`

## Benchmark Suites & Scope

1. **Univariate Summary Statistics** (`univ`, 9 datasets)
   - Datasets: `PiDigits`, `Lottery`, `Lew`, `Mavro`, `Michelso`, `NumAcc1`, `NumAcc2`, `NumAcc3`, `NumAcc4`
   - Certified statistics: Sample Mean, Sample Standard Deviation ($N-1$ denominator), Lag-1 Sample Autocorrelation ($r_1$).
2. **Linear Regression** (`lls`, 11 datasets)
   - Datasets: `Norris`, `Pontius`, `NoInt1`, `NoInt2`, `Filip`, `Longley`, `Wampler1`–`Wampler5`
   - Certified statistics: Minimum LRE across all regression coefficients (`Coef.`), minimum LRE across all standard errors (`S.E.`), and LRE of the residual sum of squares (`RSS`).
   - Note on `NoInt1` and `NoInt2`: Fit without an intercept term.
   - Note on `Filip` (Filippelli): Degree-10 polynomial with extreme multicollinearity; track both standard OLS solver behavior and orthogonal-polynomial/SVD solutions.
3. **Analysis of Variance — One-Way ANOVA** (`anova`, 11 base datasets + 3 scaled variants)
   - Base datasets: `SiRstv`, `SmLs01`–`SmLs09`, `AtmWtAg`
   - Scaled variants (`SmLs07b`, `SmLs08b`, `SmLs09b`): Multiply response $y$ by 10 before IEEE-754 binary64 storage to eliminate decimal-to-binary representation truncation (per Stata's certification methodology).
   - Certified statistic: $F$-statistic.
4. **Nonlinear Least Squares Regression** (`nls`, 27 datasets)
   - Datasets: `Misra1a`–`Misra1d`, `Chwirut1`–`2`, `Lanczos1`–`3`, `Gauss1`–`3`, `DanWood`, `Kirby2`, `Hahn1`, `Nelson`, `MGH09`, `MGH10`, `MGH17`, `Roszman1`, `ENSO`, `Thurber`, `BoxBOD`, `Rat42`, `Rat43`, `Eckerle4`, `Bennett5`
   - Certified statistics: Minimum LRE across parameters (`Coef.`), minimum LRE across standard errors (`S.E.`), and LRE of `RSS`.
   - Starting values: Attempt `Start 1` first; if convergence fails or lands on a non-optimal local minimum, fall back to `Start 2` and denote `Start 2` results in parentheses `(x.x)` just like the Stata report.

## Numerical Accuracy Metric (LRE)

Let $c$ be the computed statistic and $t$ be the NIST certified value. The **Log Relative Error (LRE)** ("Decimal Digits of Accuracy") is:

$$\text{LRE}(c, t) = \begin{cases} \min\!\left(15.0,\; -\log_{10}\!\left(\frac{|c - t|}{|t|}\right)\right) & \text{if } t \neq 0 \\ \min\!\left(15.0,\; -\log_{10}(|c - t|)\right) & \text{if } t = 0 \end{cases}$$

- Negative values are floored at `0.0` (or reported as failed convergence if applicable).
- If $c = t$ within floating-point precision, $\text{LRE} = 15.0$.

## Architecture & Directory Layout

- `data/` — Cached NIST `.dat` files and Stata reference results (`stata17_reference.json`), organized by suite (`univ/`, `lls/`, `anova/`, `nls/`).
- `src/` (or modular scripts):
  - `scraper.py` — Crawls the NIST StRD archive pages to download all `.dat` files and scrapes the Stata 17 certification page (including `.do` and `.log` links) for side-by-side comparison.
  - `parser.py` — Parses NIST StRD `.dat` files into structured metadata (difficulty, model specification, certified values, starting values, and raw observations).
  - `runners/` — Suite runners (`univariate.py`, `linear_regression.py`, `anova.py`, `nonlinear_regression.py`) implementing the statistical estimators in Python (`numpy`, `scipy`, `statsmodels`).
  - `lre.py` — LRE calculation and aggregation helpers.
  - `report.py` — Generates console summary tables, per-test logs, and an HTML/Markdown certification report mirroring the Stata 17 report alongside Stata 17 comparison numbers.
- `run_certification.py` — CLI entry point to scrape (if not cached), parse, run all tests, and generate the certification report.
