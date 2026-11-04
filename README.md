# NIST StRD Numerical Accuracy Certification in Python (vs. Stata 17)

This repository replicates the **National Institute of Standards and Technology (NIST) Statistical Reference Datasets (StRD)** numerical accuracy benchmark in Python, modeled after the [Stata 17 NIST StRD certification report](https://www.stata.com/support/cert/nist17/).

- **NIST StRD Data Archive**: <https://www.itl.nist.gov/div898/strd/general/dataarchive.html>
- **Stata 17 Certification Reference**: <https://www.stata.com/support/cert/nist17/>
- **Full Generated Reports**:
  - Markdown Report: [`reports/certification_report.md`](reports/certification_report.md)
  - Standalone HTML Report: [`reports/index.html`](reports/index.html)
  - Machine-Readable Results: [`reports/results.json`](reports/results.json)
  - Per-Test Detailed Logs: [`reports/logs/`](reports/logs/)

---

## 1. Methodology & Accuracy Metric (LRE)

NIST provides 58 benchmark datasets across four statistical suites with certified reference values computed using 500-digit multiple-precision arithmetic. Following Stata's certification methodology, numerical accuracy is measured by the **Log Relative Error (LRE)** ("Decimal Digits of Accuracy"), capped at `15.0` (IEEE-754 double precision):

$$\text{LRE}(c, t) = \begin{cases} \min\!\left(15.0,\; -\log_{10}\!\left(\frac{|c - t|}{|t|}\right)\right) & \text{if } |t| \neq 0 \\ \min\!\left(15.0,\; -\log_{10}(|c - t|)\right) & \text{otherwise} \end{cases}$$

where $c$ is the computed statistic and $t$ is the NIST certified value. For multi-parameter regression models (Linear and Nonlinear Regression), reported `Coef.` and `S.E.` values represent the **minimum LRE across all model parameters** (including the intercept, if present).

---

## 2. Executive Summary of Findings

Results were obtained running Python 3.14.4 (`numpy 2.5.4`, `scipy 1.18.1`, `statsmodels 0.15.0`, `mpmath 1.4.1`) on Linux x86\_64.

| Benchmark Suite | Certified Statistic | Python Primary (`float64`) Avg (Min–Max) | Python Extended (`50-Digit`) Avg (Min–Max) | Stata 17 (`float64`) Avg (Min–Max) | Key Finding |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **1. Univariate** (9 datasets) | Mean<br>Std. Dev. ($N-1$)<br>Lag-1 Autocorr. | **15.0** (15.0–15.0)<br>**13.3** (8.3–15.0)<br>**13.8** (10.7–15.0) | **15.0** (15.0–15.0)<br>**15.0** (15.0–15.0)<br>**15.0** (14.8–15.0) | **15.0** (15.0–15.0)<br>**13.3** (8.3–15.0)<br>**13.8** (10.7–15.0) | Two-pass `float64` matches Stata 17; naive 1-pass collapses to `0.0` on `NumAcc-4`. |
| **2. Linear Regression** (11 datasets) | Coef. (excl. Filip)<br>S.E. (excl. Filip)<br>RSS (excl. Filip)<br>`Filippelli` | **11.4** (6.7–15.0)<br>**13.4** (10.3–15.0)<br>**14.4** (12.4–15.0)<br>**7.8 / 7.8 / 9.4** | **14.9** (14.4–15.0)<br>**14.8** (14.5–15.0)<br>**15.0** (14.7–15.0)<br>**14.3 / 14.7 / 15.0** | **10.3** (6.4–15.0)<br>**13.2** (10.8–15.0)<br>**14.3** (12.7–15.0)<br>*no full solution* | Column-pivoted QR (`scipy.linalg.qr`) solves all 11 datasets including degree-10 `Filippelli` without forming $X'X$. |
| **3. One-Way ANOVA** (14 tests) | $F$-statistic (excl. unscaled S-L 7–9) | **13.1** (10.2–15.0) | **14.9** (14.7–15.0)* | **12.8** (10.2–15.0) | `scipy.stats.f_oneway` achieves `15.0` on `Simon-Lesage 2, 3, 7b, 8b, 9b`; 50-digit `Decimal` achieves `15.0` even on **unscaled** `S-L 7–9`. |
| **4. Nonlinear Regression** (27 datasets) | Coef.<br>S.E.<br>RSS | **9.3** (7.1–11.4)<br>**9.1** (3.1–10.8)<br>**10.9** (2.8–12.2) | **10.6** (7.9–12.6)**<br>**10.5** (8.3–11.2)**<br>**11.0** (10.4–12.2)** | **7.8** (4.7–10.6)<br>**5.8** (3.3–6.7)<br>**10.9** (3.0–12.2) | Converges from `Start 1` on all 27 datasets (Stata required `Start 2` on 5 datasets) and gains +3.3 digits on S.E. via complex-step Jacobians. |

*\*In 50-digit `Decimal`, unscaled `Simon-Lesage 7`, `8`, and `9` also achieve **15.0** digits of accuracy.*  
*\*\*NIST `.dat` files for Nonlinear Regression publish certified values rounded to 11 significant digits ($\approx 10.5\text{–}11.2$ LRE ceiling).*

---

## 3. Detailed Certification Results by Suite

### 3.1 Univariate Summary Statistics

#### Primary Run (`float64` `numpy` / `statsmodels`) vs. Stata 17

```text
                                          Python                   Stata 17
                                    Digits of accuracy        Digits of accuracy
                                  ----------------------    ----------------------
                                                  lag-1                     lag-1
        Test           Difficulty  mean    S.D.  autocorr.   mean    S.D.  autocorr.
        ----------------------------------------------------------------------------
        PiDigits       lower       15.0    15.0      14.0    15.0    15.0      14.9
        Lottery        lower       15.0    15.0      14.9    15.0    15.0      15.0
        Lew            lower       15.0    15.0      14.8    15.0    15.0      14.8
        Mavro          lower       15.0    13.1      14.1    15.0    13.1      13.7
        Michelson      lower       15.0    13.8      13.4    15.0    13.8      13.4
        NumAcc-1       lower       15.0    15.0      15.0    15.0    15.0      15.0
        NumAcc-2       average     15.0    15.0      14.9    15.0    15.0      15.0
        NumAcc-3       average     15.0     9.5      11.9    15.0     9.5      11.9
        NumAcc-4       higher      15.0     8.3      10.7    15.0     8.3      10.7
        ----------------------------------------------------------------------------
        Average                    15.0    13.3      13.8    15.0    13.3      13.8
        Minimum                    15.0     8.3      10.7    15.0     8.3      10.7
        Maximum                    15.0    15.0      15.0    15.0    15.0      15.0
```

#### Additional Methods: Algorithm & Extended-Precision Comparison

```text
                       Primary (Two-Pass)    1-Pass Textbook    80-bit longdouble    50-Digit Decimal       Stata 17
                            (float64)           (float64)         (np.float128)      (Exact Decimal)       (float64)
                       ------------------   ----------------   ------------------   ------------------   --------------
        Test             S.D.   autocorr.    S.D.  autocorr.     S.D.   autocorr.     S.D.   autocorr.    S.D. autocorr.
        ---------------------------------------------------------------------------------------------------------------
        PiDigits         15.0       14.0     15.0      15.0      15.0       15.0      15.0       15.0     15.0     14.9
        Lottery          15.0       14.9     15.0      15.0      15.0       14.9      15.0       14.9     15.0     15.0
        Lew              15.0       14.8     15.0      14.8      15.0       14.8      15.0       14.8     15.0     14.8
        Mavro            13.1       14.1      8.2       7.9      15.0       15.0      15.0       15.0     13.1     13.7
        Michelson        13.8       13.4      8.6       8.3      15.0       15.0      15.0       15.0     13.8     13.4
        NumAcc-1         15.0       15.0     15.0      15.0      15.0       15.0      15.0       15.0     15.0     15.0
        NumAcc-2         15.0       14.9     13.6      13.3      15.0       15.0      15.0       15.0     15.0     15.0
        NumAcc-3          9.5       11.9      1.9       1.6      13.2       15.0      15.0       15.0      9.5     11.9
        NumAcc-4          8.3       10.7      0.0       0.0      12.0       15.0      15.0       15.0      8.3     10.7
        ---------------------------------------------------------------------------------------------------------------
        Average          13.3       13.8     10.3      10.1      14.5       15.0      15.0       15.0     13.3     13.8
```

---

### 3.2 Linear Regression

#### Primary Run (`scipy.linalg` Column-Pivoted QR, `float64`) vs. Stata 17

```text
                                          Python                   Stata 17
                                    Digits of accuracy        Digits of accuracy
                                  ----------------------    ----------------------
        Test           Difficulty  Coef.   S.E.      RSS     Coef.   S.E.      RSS
        --------------------------------------------------------------------------
        Norris         lower       13.2    13.9     13.7     12.8    13.5     13.3
        Pontius        lower       12.3    13.6     13.3     11.5    13.0     12.7
        NoInt-1        average     14.7    15.0     14.9     14.7    15.0     14.9
        NoInt-2        average     15.0    14.9     15.0     15.0    15.0     14.7
        Filippelli     higher       7.8     7.8      9.4      no full solution(*)
        Longley        higher      11.1    12.7     12.4     12.1    12.9     13.2
        Wampler-1      higher       9.9    10.3     15.0      6.9    15.0     15.0
        Wampler-2      higher      12.7    14.4     15.0     10.4    15.0     15.0
        Wampler-3      higher      10.1    13.2     14.2      6.5    10.8     14.1
        Wampler-4      higher       8.7    13.1     15.0      6.5    10.8     15.0
        Wampler-5      higher       6.7    13.1     15.0      6.4    10.8     15.0
        --------------------------------------------------------------------------
        Average (all 11)           11.1    12.9     13.9        -       -        -
        Average (excl. Filip)      11.4    13.4     14.4     10.3    13.2     14.3
        Minimum (excl. Filip)       6.7    10.3     12.4      6.4    10.8     12.7
        Maximum                    15.0    15.0     15.0     15.0    15.0     15.0
        --------------------------------------------------------------------------
```

#### Additional Methods: Multi-Solver & Extended-Precision Comparison

- **Why Column-Pivoted QR + Triangular Solve wins in `float64`**:
  - Forming normal equations $(X'X)^{-1}X'y$ squares the condition number $\kappa(X)^2$, causing `Filippelli` ($\kappa(X) \approx 1.2 \times 10^{15}$) to fail completely and degrading `Wampler-1..5`.
  - `statsmodels.OLS(method="pinv")` (default SVD) truncates small singular values below machine-epsilon threshold, marking `Filippelli` as rank-deficient (`n/a`, similar to Stata's `regress`) and losing accuracy on `Pontius` (`6.2`).
  - `statsmodels.OLS(method="qr")` computes coefficients via QR (`8.5` on `Filippelli`), but computes standard errors via `np.linalg.inv(R.T @ R)`, which squares $\kappa(R)$ and produces negative diagonal entries (`NaN` S.E.) on `Filippelli`.
  - Inverting $R$ directly via `scipy.linalg.solve_triangular(R, I)` and summing row norms (`np.sum(R_inv**2, axis=1)`) avoids squaring $\kappa(R)$ and yields `7.8` digits on both coefficients and standard errors for `Filippelli`.

```text
                       Primary: Pivoted QR   statsmodels QR    statsmodels SVD   Normal Equations   50-Digit mpmath      Stata 17
                         (scipy float64)        (float64)      (pinv, float64)    (X'X, float64)    (Exact 50-dps)      (float64)
                       -------------------   --------------   ----------------   ----------------   ---------------   --------------
        Test              Coef.     S.E.      Coef.    S.E.     Coef.     S.E.     Coef.     S.E.     Coef.    S.E.    Coef.    S.E.
        ----------------------------------------------------------------------------------------------------------------------------
        Norris             13.2     13.9       13.3    13.8      12.8     13.6      12.1     13.8      14.4    14.7     12.8    13.5
        Pontius            12.3     13.6       12.4    13.8       6.2      9.2      11.2     13.7      15.0    14.7     11.5    13.0
        NoInt-1            14.7     15.0       14.7    15.0      14.7     15.0      14.7     15.0      14.7    15.0     14.7    15.0
        NoInt-2            15.0     14.9       15.0    14.9      15.0     14.8      15.0     15.0      15.0    14.9     15.0    15.0
        Filippelli          7.8      7.8        8.5     n/a       n/a      n/a       0.0      n/a      14.3    14.7      n/a     n/a
        Longley            11.1     12.7       10.9     8.0      11.0     11.8       7.4      8.5      14.6    14.8     12.1    12.9
        Wampler-1           9.9     10.3        9.5     9.7       9.6      8.3       6.2      6.8      15.0    15.0      6.9    15.0
        Wampler-2          12.7     14.4       13.2    14.3      10.1     10.6       9.4     11.6      15.0    15.0     10.4    15.0
        Wampler-3          10.1     13.2        9.2    10.4       9.4      9.8       6.2     10.5      15.0    14.5      6.5    10.8
        Wampler-4           8.7     13.1        7.7    10.4       7.7      9.8       6.2     10.5      15.0    14.5      6.5    10.8
        Wampler-5           6.7     13.1        5.7    10.4       5.7      9.8       6.2     10.5      15.0    14.5      6.4    10.8
        ----------------------------------------------------------------------------------------------------------------------------
```

---

### 3.3 Analysis of Variance (One-Way ANOVA)

#### Primary Run (`scipy.stats.f_oneway`, `float64`) vs. Stata 17

```text
                                              Python             Stata 17
                                        Digits of accuracy  Digits of accuracy
                                        ------------------  ------------------
        Test               Difficulty            F                   F
        ----------------------------------------------------------------------
        Si Resistivity     lower                  13.1                13.1
        Simon-Lesage 1     lower                  15.0                15.0
        Simon-Lesage 2     lower                  15.0                13.6
        Simon-Lesage 3     lower                  15.0                12.8
        Ag Atomic Wt       average                10.2                10.2
        Simon-Lesage 4     average                10.4                10.4
        Simon-Lesage 5     average                10.2                10.2
        Simon-Lesage 6     average                10.2                10.2
        Simon-Lesage 7     higher               4.4(*)              4.4(*)
                     7b    higher              15.0(*)             15.0(*)
        Simon-Lesage 8     higher               4.2(*)              4.3(*)
                     8b    higher              15.0(*)             15.0(*)
        Simon-Lesage 9     higher               4.2(*)              4.2(*)
                     9b    higher              15.0(*)             15.0(*)
        ----------------------------------------------------------------------
        Average excluding S-L 7, 8, 9          13.1                12.8
        Minimum                                10.2                10.2
        Maximum                                15.0                15.0
        ----------------------------------------------------------------------
```

#### Additional Methods: Multi-Solver & Extended-Precision Comparison

- **Confirming the Simon–Lesage 7–9 Binary64 Truncation Effect**:
  - In `Simon-Lesage 7`, `8`, and `9`, observations take values such as `1000000000000.4`, which cannot be represented exactly in 53-bit binary64 (`1000000000000.4000244140625`).
  - Multiplying inputs by 10 in decimal before binary64 conversion (`7b`, `8b`, `9b`) produces exact integers (`10000000000004`), allowing `scipy.stats.f_oneway` to achieve **15.0** digits of accuracy.
  - Running the **unscaled** `Simon-Lesage 7`, `8`, and `9` datasets in 80-bit `np.longdouble` gains ~3.3 digits (`7.5`–`7.7`), and running them in 50-digit `Decimal` achieves **15.0** digits without scaling.

```text
                            Primary:       2-Pass     statsmodels    1-Pass       80-bit      50-Digit     Stata 17
                          scipy.stats     Centered        OLS       Textbook    longdouble    Decimal     Reference
                           (float64)      (float64)    (float64)    (float64)   (float128)   (50-Digit)   (float64)
        Test                   F              F            F            F            F            F           F
        -----------------------------------------------------------------------------------------------------------
        Si Resistivity        13.1           13.1         13.1          9.5         14.7         14.7        13.1
        Simon-Lesage 1        15.0           15.0         15.0         12.6         15.0         15.0        15.0
        Simon-Lesage 2        15.0           15.0         15.0         13.3         15.0         15.0        13.6
        Simon-Lesage 3        15.0           15.0         14.6         12.9         15.0         15.0        12.8
        Ag Atomic Wt          10.2           10.2         10.2          1.4         13.0         14.7        10.2
        Simon-Lesage 4        10.4           10.4         10.4          1.8         13.7         15.0        10.4
        Simon-Lesage 5        10.2           10.2         10.2          1.1         13.5         15.0        10.2
        Simon-Lesage 6        10.2           10.2         10.2          2.7         13.5         15.0        10.2
        Simon-Lesage 7         4.4            4.4          5.4          0.0          7.7         15.0         4.4
        Simon-Lesage 7b       15.0           15.0          4.0          0.0         15.0         15.0        15.0
        Simon-Lesage 8         4.2            4.2          2.4          0.9          7.5         15.0         4.3
        Simon-Lesage 8b       15.0           15.0          1.9          0.0         15.0         15.0        15.0
        Simon-Lesage 9         4.2            4.2          1.9          0.0          7.5         15.0         4.2
        Simon-Lesage 9b       15.0            5.4          1.6          0.0         15.0         15.0        15.0
        -----------------------------------------------------------------------------------------------------------
```

---

### 3.4 Nonlinear Least Squares Regression

#### Primary Run (`scipy.optimize.least_squares` LM + Complex-Step Jacobian, `float64`) vs. Stata 17

```text
                                           Python                      Stata 17
                                     Digits of accuracy           Digits of accuracy
                                  -------------------------    -------------------------
        Test           Difficulty   Coef.    S.E.      RSS       Coef.    S.E.      RSS
        --------------------------------------------------------------------------------
        Misra 1a       lower          9.8     9.5     10.5         9.4     6.4     10.5
        Chwirut 2      lower         10.1    10.8     11.2         8.0     6.3     11.2
        Chwirut 1      lower          9.7    10.0     11.4         7.6     6.3     11.4
        Lanczos 3      lower          8.2     8.2     10.6         7.2     6.0     10.6
        Gauss 1        lower         10.6    10.6     11.6         8.5     6.3     11.6
        Gauss 2        lower         10.0     9.9     10.6         8.2     5.9     10.6
        Daniel Wood    lower         10.0    10.1     11.7         8.6     6.2     11.7
        Misra 1b       lower          9.5     9.3     11.3         9.9     6.5     11.3
        Kirby 2        average        7.7     8.1     11.6         8.0     6.3     11.6
        Hahn 1         average        7.5     8.3     10.6         7.1     5.1     10.6
        Nelson         average        7.9     8.0     10.9         7.1     5.2     10.9
        MGH 17         average        8.1     7.7     11.5       (7.0)   (6.1)   (11.5)
        Lanczos 1      average       10.6     3.1      2.8        10.6     3.3      3.0
        Lanczos 2      average       10.4    10.2     10.0         7.9     5.4     10.1
        Gauss 3        average       10.4    10.3     11.0         8.2     5.5     11.0
        Misra 1c       average       10.8    10.6     11.1         9.7     6.5     11.1
        Misra 1d       average       11.4    10.8     11.2         9.3     6.5     11.2
        Roszman 1      average        7.1     7.7     12.2         7.4     6.4     12.2
        ENSO           average        7.5     8.2     11.4         4.7     5.3     11.3
        MGH 09         higher         7.2     7.3     11.6       (7.0)   (6.5)   (11.6)
        Thurber        higher         8.2     7.7     11.2         6.5     5.4     11.3
        BoxBOD         higher         9.3     9.3     10.4         7.3     6.7     10.4
        Ratkowsky 2    higher        10.4    10.1     11.9         7.6     6.0     11.8
        MGH 10         higher         9.9    10.1     11.4       (7.7)   (4.7)   (11.4)
        Eckerle4       higher         9.7     9.7     10.7       (8.3)   (6.4)   (10.7)
        Ratkowsky 3    higher         8.0     8.2     11.4       (6.0)   (5.0)   (11.4)
        Bennett 5      higher        11.0    10.4     11.2         6.4     5.9     11.0
        --------------------------------------------------------------------------------
        Average                       9.3     9.1     10.9         7.8     5.8     10.9
        Minimum                       7.1     3.1      2.8         4.7     3.3      3.0
        Maximum                      11.4    10.8     12.2        10.6     6.7     12.2
        --------------------------------------------------------------------------------
```

#### Additional Methods: Multi-Solver & Extended-Precision Comparison

- **Complex-Step Differentiation vs. Finite-Difference Jacobians**:
  - When `scipy.optimize.least_squares` is run with default 2-point finite differences (`LM + 2-pt Finite Diff`), it averages **7.5** digits on `Coef.` and **5.5** digits on `S.E.`, closely mirroring Stata 17's **7.8** `Coef.` and **5.8** `S.E.`. This demonstrates that Stata's ~5.8-digit S.E. ceiling is driven by numerical derivative step-size truncation ($\sim \sqrt{\epsilon_{\text{mach}}}$).
  - Evaluating the Jacobian via complex-step differentiation $\text{Im}(f(b + ih e_j))/h$ eliminates subtractive cancellation, raising `float64` S.E. accuracy from **5.5** to **9.1** digits and achieving convergence from `Start 1` on all 27 datasets.
  - On `Lanczos 1` (exact synthetic zero-residual data where true $\text{RSS} \approx 1.43 \times 10^{-25}$), `float64` residuals hit the $10^{-16}$ rounding floor (`2.8`–`3.3` digits in both Python `float64` and Stata 17), whereas 50-digit `mpmath` recovers the full **10.6** digits of the certified table.

```text
                       Primary: LM + Complex-Step    LM + 2-pt Finite Diff    TRF + 3-pt Central Diff    50-Digit mpmath         Stata 17
                           (scipy, float64)            (scipy, float64)          (scipy, float64)        (Exact 50-dps)          (float64)
                       --------------------------    ---------------------    -----------------------    ---------------    ------------------
        Test                Coef.        S.E.          Coef.        S.E.          Coef.        S.E.       Coef.     S.E.      Coef.      S.E.
        --------------------------------------------------------------------------------------------------------------------------------------
        Misra 1a               9.8         9.5            7.4         4.6            6.9         4.7       11.1     10.8        9.4       6.4
        Chwirut 2             10.1        10.8            9.1         5.6            8.7         6.2       10.9     11.0        8.0       6.3
        Chwirut 1              9.7        10.0            8.3         5.6            8.4         6.3       10.7     10.9        7.6       6.3
        Lanczos 3              8.2         8.2            5.7         4.3            8.9         7.7       10.5     11.2        7.2       6.0
        Gauss 1               10.6        10.6            8.1         5.9            8.3         6.8       10.6     10.6        8.5       6.3
        Gauss 2               10.0         9.9            9.2         5.9            9.6         6.7       10.3     10.4        8.2       5.9
        Daniel Wood           10.0        10.1            9.3         7.8           10.2        10.1       11.2     11.1        8.6       6.2
        Misra 1b               9.5         9.3            7.3         4.5            6.6         4.3       11.0     10.8        9.9       6.5
        Kirby 2                7.7         8.1            5.1         3.1            2.5         1.1       10.5     10.7        8.0       6.3
        Hahn 1                 7.5         8.3            2.2         0.4          (0.2)       (0.0)       10.8     10.5        7.1       5.1
        Nelson                 7.9         8.0            7.4         5.4            7.4         5.9       10.9     10.6        7.1       5.2
        MGH 17                 8.1         7.7            6.4         5.4            7.6         6.0       10.8     10.6      (7.0)     (6.1)
        Lanczos 1             10.6         3.1           10.6         3.0           10.6         3.3       10.6     10.6       10.6       3.3
        Lanczos 2             10.4        10.2            7.1         4.5            9.2         7.4       10.4     10.8        7.9       5.4
        Gauss 3               10.4        10.3            8.5         5.9            9.3         6.7       10.5     10.5        8.2       5.5
        Misra 1c              10.8        10.6            7.1         4.2            6.3         3.7       10.8     10.6        9.7       6.5
        Misra 1d              11.4        10.8            7.2         4.4            6.5         4.0       12.6     11.1        9.3       6.5
        Roszman 1              7.1         7.7            8.3         7.0            9.0         9.6       10.9     11.0        7.4       6.4
        ENSO                   7.5         8.2            6.1         6.7            6.6         7.7        7.9      9.1        4.7       5.3
        MGH 09                 7.2         7.3            7.4         7.0            7.4         7.4        8.2      8.3      (7.0)     (6.5)
        Thurber                8.2         7.7            7.2         6.4            7.7         6.8        8.9      8.3        6.5       5.4
        BoxBOD                 9.3         9.3          (8.3)       (7.6)            9.0         9.0       11.7     10.4        7.3       6.7
        Ratkowsky 2           10.4        10.1            7.8         6.5           10.8         7.6       11.2     10.6        7.6       6.0
        MGH 10                 9.9        10.1            7.5         6.3            9.4         8.0       10.9     11.1      (7.7)     (4.7)
        Eckerle4               9.7         9.7            9.6         7.8            9.3         6.9       10.5     10.5      (8.3)     (6.4)
        Ratkowsky 3            8.0         8.2            7.4         6.4            8.5         8.6       11.0     10.7      (6.0)     (5.0)
        Bennett 5             11.0        10.4            5.8         5.0            7.7         7.7       11.0     10.4        6.4       5.9
        --------------------------------------------------------------------------------------------------------------------------------------
        Average                9.3         9.1            7.5         5.5            7.9         6.3       10.6     10.5        7.8       5.8
        Minimum                7.1         3.1            2.2         0.4            0.2         0.0        7.9      8.3        4.7       3.3
        --------------------------------------------------------------------------------------------------------------------------------------
```

---

## 4. Replication Instructions

### 4.1 Environment Setup

Create a virtual environment and install the dependencies listed in [`requirements.txt`](requirements.txt):

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

### 4.2 Run the Scraper & Certification Benchmark

Run [`run_certification.py`](run_certification.py) from the repository root:

```bash
# Uses cached NIST .dat files and Stata 17 reference data in ./data (or scrapes them if not present)
.venv/bin/python run_certification.py

# Force a fresh scrape of all NIST StRD web pages and Stata 17 certification files
.venv/bin/python run_certification.py --refresh
```

This command:
1. Crawls the [NIST StRD archive](https://www.itl.nist.gov/div898/strd/general/dataarchive.html) to download all 58 `.dat` files across `univ`, `lls`, `anova`, and `nls` into `data/nist/`.
2. Scrapes the [Stata 17 certification page](https://www.stata.com/support/cert/nist17/) and downloads all 61 Stata `.do` files into `data/stata17/`.
3. Executes the primary `float64` statistical runners alongside the multi-solver and extended-precision (`np.longdouble`, `Decimal`, and `mpmath` 50-digit) runners.
4. Prints the full certification tables to the terminal and writes:
   - `reports/certification_report.md`
   - `reports/index.html`
   - `reports/results.json`
   - `reports/logs/<suite>/<dataset>.log` (61 individual test log files)

### 4.3 Run the Automated Test Suite

Verify the LRE evaluator, `.dat` parser, primary `float64` runners, and 50-digit extended-precision runners via `pytest`:

```bash
.venv/bin/python -m pytest -v
```
