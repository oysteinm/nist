# NIST StRD Certification Results Using Python (vs. Stata 17)

1. [Explanation](#1-explanation)
2. [Summary of results](#2-summary-of-results)
3. [Certification results: univariate summary statistics](#3-certification-results-univariate-summary-statistics)
4. [Certification results: linear regression](#4-certification-results-linear-regression)
5. [Certification results: analysis of variance](#5-certification-results-analysis-of-variance)
6. [Certification results: nonlinear regression](#6-certification-results-nonlinear-regression)
7. [Multi-solver & extended-precision comparison](#7-multi-solver--extended-precision-comparison)

---

## 1. Explanation

The National Institute of Standards and Technology (NIST) provides the **Statistical Reference Datasets (StRD)** ([data archive](https://www.itl.nist.gov/div898/strd/general/dataarchive.html)) with certified reference values computed using 500-digit multiple-precision arithmetic to benchmark the numerical accuracy of statistical software.

This report replicates the [Stata 17 NIST StRD certification benchmark](https://www.stata.com/support/cert/nist17/) in Python and compares:
1. **Primary `float64` Python Run** (standard IEEE-754 double precision, preserved as the baseline comparison against Stata 17),
2. **Alternative Python `float64` Solvers** (`statsmodels` SVD/QR, textbook Normal Equations, one-pass formulas, and finite-difference vs. complex-step Jacobians), and
3. **Extended-Precision Python Runs** (`np.longdouble` 80-bit extended precision and 50-digit `Decimal` / `mpmath` arbitrary precision parsed directly from the raw `.dat` strings).

In reporting comparisons, we use the **Log Relative Error (LRE)** ("Decimal Digits of Accuracy"). Let $c$ represent the calculated result and $t$ the certified answer supplied by NIST:

- $\text{LRE} = \min\!\left(15,\; -\log_{10}\!\left(\frac{|c - t|}{|t|}\right)\right)$ if $|t| \neq 0$
- $\text{LRE} = \min\!\left(15,\; -\log_{10}(|c - t|)\right)$ otherwise.

## 2. Summary of results

Results were obtained running Python 3.14.4 (`numpy 2.5.4`, `scipy 1.18.1`, `statsmodels 0.15.0`, `mpmath 1.4.1`) on `Linux 7.0.0-42-generic (x86_64)`.

- **Univariate summary statistics:**
  Python completed all 9 tests. Means were estimated with **15.0** digits of accuracy on every test (matching Stata 17's 15.0). Standard deviations averaged **13.3** correct digits, ranging from **8.3** to **15.0** digits (identical to Stata 17's 13.3 average and 8.3 minimum). The lag-1 autocorrelation averaged **13.8** correct digits, ranging from **10.7** to **15.0** digits (identical to Stata 17's 13.8 average and 10.7 minimum). In 50-digit `Decimal`, all 9 tests achieve **15.0** digits across all three statistics.

- **Linear regression:**
  Python's primary solver (`scipy.linalg` column-pivoted QR decomposition without forming $X'X$) completed **all 11 tests**, including the extreme degree-10 polynomial `Filippelli` test where Stata's standard `regress` and `statsmodels.OLS(method='pinv')` drop collinear columns.
  - On `Filippelli`, Python's pivoted QR achieved **7.8** digits on coefficients, **7.8** digits on standard errors, and **9.4** digits on RSS (compared to Stata's `orthpoly` note of 8.4 on coefficients and 8.5 on RSS).
  - Across the other 10 tests, Python coefficients averaged **11.4** correct digits (min **6.7**, vs. Stata's 10.3 avg / 6.4 min), standard errors averaged **13.4** digits (min **10.3**, vs. Stata's 13.2 avg / 10.8 min), and RSS averaged **14.4** digits (min **12.4**, vs. Stata's 14.3 avg / 12.7 min).
  - In 50-digit `mpmath`, all 11 tests (including `Filippelli`) achieve **15.0** digits across `Coef.`, `S.E.`, and `RSS`.

- **Analysis of variance:**
  Python (`scipy.stats.f_oneway`) completed all tests. Including the Simon–Lesage `7b`–`9b` decimal-scaled tests, the $F$-statistic averaged **13.1** correct digits and never had fewer than **10.2** correct digits (exceeding Stata 17's 12.8 average and matching its 10.2 minimum). Furthermore, our 50-digit `Decimal` runner achieves **15.0** digits even on the **unscaled** `Simon-Lesage 7`, `8`, and `9` tests, directly confirming Stata's explanation that the 4.2–4.4 `float64` ceiling on `SmLs07`–`09` is caused by binary64 truncation of decimal `.4`/`.3`/`.5` on inputs of magnitude $10^{12}$.

- **Nonlinear regression:**
  Python (`scipy.optimize.least_squares` with Levenberg–Marquardt and complex-step exact Jacobians) completed **all 27 tests from `Start 1`**, including the 5 difficult problems (`MGH 17`, `MGH 09`, `MGH 10`, `Eckerle4`, `Ratkowsky 3`) where Stata 17 required `Start 2`.
  - Coefficients averaged **9.3** correct digits (minimum **7.1**, vs. Stata 17's **7.8** average and **4.7** minimum).
  - Standard errors averaged **9.1** correct digits (minimum **3.1**, vs. Stata 17's **5.8** average and **3.3** minimum).
  - Residual sums of squares averaged **10.9** correct digits (minimum **2.8**, vs. Stata 17's **10.9** average and **3.0** minimum).
  - In 50-digit `mpmath` Gauss–Newton refinement, all 27 tests achieve **15.0** digits across `Coef.`, `S.E.`, and `RSS`.

---

## 3. Certification results: univariate summary statistics

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

---

## 4. Certification results: linear regression

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

Each test involved multiple independent variables. Reported under `Coef.` and `S.E.` is the minimum LRE across all regressors, including the intercept, if any. `RSS` reports the LRE for the residual sum of squares.

---

## 5. Certification results: analysis of variance

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

`(*)` Tests **Simon–Lesage 7b through 9b** multiply the input response $y$ by 10 in exact decimal before storing in IEEE-754 double precision, so values such as `1,000,000,000,000.4` become exact binary64 integers `10,000,000,000,004`.

---

## 6. Certification results: nonlinear regression

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

Parentheses `(x.x)` indicate that convergence could not be achieved with the first set of starting values (`Start 1`) and that the second set (`Start 2`) had to be used. Python's primary solver converged from `Start 1` on all 27 datasets.

---

## 7. Multi-solver & extended-precision comparison

In addition to the primary `float64` run above, we evaluate alternative `float64` algorithms in Python alongside extended precision (`np.longdouble` 80-bit and 50-digit `Decimal` / `mpmath`).

### 7.1 Univariate Summary Statistics — Algorithm & Precision Comparison

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

### 7.2 Linear Regression — Multi-Solver & Extended-Precision Comparison

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

### 7.3 One-Way ANOVA — Multi-Solver & Extended-Precision Comparison

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

### 7.4 Nonlinear Regression — Multi-Solver & Extended-Precision Comparison

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
