"""Linear Least Squares Regression runner for NIST StRD datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import warnings

import mpmath as mp
import numpy as np
import scipy.linalg
import statsmodels.api as sm

from src.lre import compute_lre, min_lre
from src.parser import LinearRegressionDataset, parse_lls_dat

LLS_ORDER = [
    ("Norris", "Norris"),
    ("Pontius", "Pontius"),
    ("NoInt1", "NoInt-1"),
    ("NoInt2", "NoInt-2"),
    ("Filip", "Filippelli"),
    ("Longley", "Longley"),
    ("Wampler1", "Wampler-1"),
    ("Wampler2", "Wampler-2"),
    ("Wampler3", "Wampler-3"),
    ("Wampler4", "Wampler-4"),
    ("Wampler5", "Wampler-5"),
]


@dataclass
class LinearRegressionResult:
    dataset_name: str
    display_name: str
    difficulty: str
    n_obs: int
    n_params: int
    has_intercept: bool
    param_names: list[str]
    calc_coef: list[float]
    cert_coef: list[float]
    param_lre_coef: list[float]
    lre_coef: float
    calc_se: list[float]
    cert_se: list[float]
    param_lre_se: list[float]
    lre_se: float
    calc_rss: float
    cert_rss: float
    lre_rss: float
    cond_number: float
    rank_deficient: bool = False
    note: str | None = None
    alt_solutions: dict[str, dict[str, float | None]] = field(default_factory=dict)


def build_design_matrix(ds: LinearRegressionDataset) -> np.ndarray:
    """Construct the design matrix X matching the NIST model specification."""
    n_non_intercept = len(ds.param_names) - (1 if ds.has_intercept else 0)
    if ds.X_raw.shape[1] == 1 and n_non_intercept > 1:
        x = ds.X_raw[:, 0]
        X_pred = np.column_stack([x**p for p in range(1, n_non_intercept + 1)])
    else:
        X_pred = ds.X_raw.copy()

    if ds.has_intercept:
        return np.column_stack([np.ones(len(ds.y), dtype=np.float64), X_pred])
    return X_pred


def _run_normal_equations(X: np.ndarray, y: np.ndarray, ds: LinearRegressionDataset) -> dict[str, float | None]:
    """Solve OLS via textbook Normal Equations (X'X)b = X'y in float64."""
    n, p = X.shape
    df_r = n - p
    try:
        xtx = X.T @ X
        xty = X.T @ y
        coef = np.linalg.solve(xtx, xty)
        res = y - X @ coef
        rss = float(np.dot(res, res))
        xtx_inv = np.linalg.inv(xtx)
        diag_v = np.diag(xtx_inv) * (rss / df_r)
        if np.any(diag_v < 0):
            se = None
        else:
            se = [float(v) for v in np.sqrt(diag_v)]
        return {
            "lre_coef": min_lre([float(v) for v in coef], ds.certified_coef),
            "lre_se": None if se is None else min_lre(se, ds.certified_se),
            "lre_rss": compute_lre(rss, ds.certified_rss),
        }
    except Exception:
        return {"lre_coef": None, "lre_se": None, "lre_rss": None}


def _run_mpmath_ols(ds: LinearRegressionDataset, dps: int = 50) -> dict[str, float | None]:
    """Solve OLS in 50-digit arbitrary precision via mpmath from raw .dat strings."""
    with mp.workdps(dps):
        n = len(ds.raw_rows)
        p = len(ds.param_names)
        df_r = n - p
        n_non_intercept = p - (1 if ds.has_intercept else 0)

        y_mp = mp.matrix(n, 1)
        X_mp = mp.matrix(n, p)
        for i, row in enumerate(ds.raw_rows):
            y_mp[i, 0] = mp.mpf(row[0])
            col_idx = 0
            if ds.has_intercept:
                X_mp[i, 0] = mp.mpf(1)
                col_idx = 1
            if len(row) - 1 == 1 and n_non_intercept > 1:
                x_val = mp.mpf(row[1])
                for deg in range(1, n_non_intercept + 1):
                    X_mp[i, col_idx + deg - 1] = x_val**deg
            else:
                for j in range(n_non_intercept):
                    X_mp[i, col_idx + j] = mp.mpf(row[1 + j])

        xtx = X_mp.T * X_mp
        xty = X_mp.T * y_mp
        coef_mp = mp.lu_solve(xtx, xty)
        res_mp = y_mp - X_mp * coef_mp
        rss_mp = sum(res_mp[i, 0] ** 2 for i in range(n))
        xtx_inv = xtx**-1
        se_mp = [mp.sqrt(xtx_inv[j, j] * (rss_mp / df_r)) for j in range(p)]

        calc_coef = [float(coef_mp[j, 0]) for j in range(p)]
        calc_se = [float(se_mp[j]) for j in range(p)]
        calc_rss = float(rss_mp)

        return {
            "lre_coef": min_lre(calc_coef, ds.certified_coef) or 0.0,
            "lre_se": min_lre(calc_se, ds.certified_se) or 0.0,
            "lre_rss": compute_lre(calc_rss, ds.certified_rss) or 0.0,
        }


def run_lls_dataset(ds: LinearRegressionDataset, display_name: str | None = None) -> LinearRegressionResult:
    """Fit OLS via column-pivoted QR decomposition (and alternative solvers) and compute LRE."""
    X = build_design_matrix(ds)
    y = ds.y
    n, p = X.shape
    df_r = n - p

    # Primary solver (Initial run preserved): Column-pivoted QR via scipy.linalg.qr + triangular solve
    Q, R, P = scipy.linalg.qr(X, mode="economic", pivoting=True)
    coef_piv = scipy.linalg.solve_triangular(R, Q.T @ y)
    coef_arr = np.empty_like(coef_piv)
    coef_arr[P] = coef_piv

    residuals = y - X @ coef_arr
    calc_rss = float(np.dot(residuals, residuals))

    R_inv = scipy.linalg.solve_triangular(R, np.eye(p))
    var_piv = np.sum(R_inv**2, axis=1) * (calc_rss / df_r)
    se_arr = np.empty_like(var_piv)
    se_arr[P] = np.sqrt(var_piv)

    calc_coef = [float(v) for v in coef_arr]
    calc_se = [float(v) for v in se_arr]
    cond_number = float(np.linalg.cond(X))

    param_lre_coef = [compute_lre(c, t) or 0.0 for c, t in zip(calc_coef, ds.certified_coef)]
    param_lre_se = [compute_lre(c, t) or 0.0 for c, t in zip(calc_se, ds.certified_se)]

    lre_coef = min_lre(calc_coef, ds.certified_coef) or 0.0
    lre_se = min_lre(calc_se, ds.certified_se) or 0.0
    lre_rss = compute_lre(calc_rss, ds.certified_rss) or 0.0

    # Compare statsmodels OLS (default SVD/pinv and unpivoted QR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = sm.OLS(y, X)
        res_qr = model.fit(method="qr")
        res_pinv = model.fit(method="pinv")
        rank_deficient = int(model.rank) < p
        qr_params = [float(v) for v in res_qr.params]
        qr_bse = [float(v) for v in res_qr.bse]
        qr_ssr = float(res_qr.ssr)
        pinv_params = [float(v) for v in res_pinv.params]
        pinv_bse = [float(v) for v in res_pinv.bse]
        pinv_ssr = float(res_pinv.ssr)

    alt_solutions: dict[str, dict[str, float | None]] = {
        "statsmodels_qr": {
            "lre_coef": min_lre(qr_params, ds.certified_coef),
            "lre_se": min_lre(qr_bse, ds.certified_se),
            "lre_rss": compute_lre(qr_ssr, ds.certified_rss),
        },
        "statsmodels_pinv": {
            "lre_coef": None if rank_deficient else min_lre(pinv_params, ds.certified_coef),
            "lre_se": None if rank_deficient else min_lre(pinv_bse, ds.certified_se),
            "lre_rss": None if rank_deficient else compute_lre(pinv_ssr, ds.certified_rss),
        },
        "normal_equations": _run_normal_equations(X, y, ds),
        "mpmath_50digit": _run_mpmath_ols(ds, dps=50),
    }

    note = "(*)" if ds.name == "Filip" else None

    return LinearRegressionResult(
        dataset_name=ds.name,
        display_name=display_name or ds.name,
        difficulty=ds.difficulty,
        n_obs=n,
        n_params=p,
        has_intercept=ds.has_intercept,
        param_names=ds.param_names,
        calc_coef=calc_coef,
        cert_coef=ds.certified_coef,
        param_lre_coef=param_lre_coef,
        lre_coef=lre_coef,
        calc_se=calc_se,
        cert_se=ds.certified_se,
        param_lre_se=param_lre_se,
        lre_se=lre_se,
        calc_rss=calc_rss,
        cert_rss=ds.certified_rss,
        lre_rss=lre_rss,
        cond_number=cond_number,
        rank_deficient=rank_deficient,
        note=note,
        alt_solutions=alt_solutions,
    )


def run_lls_suite(nist_lls_dir: Path) -> list[LinearRegressionResult]:
    """Run all 11 NIST StRD Linear Regression datasets in Stata report order."""
    results: list[LinearRegressionResult] = []
    for file_stem, display_name in LLS_ORDER:
        dat_path = nist_lls_dir / f"{file_stem}.dat"
        ds = parse_lls_dat(dat_path)
        results.append(run_lls_dataset(ds, display_name=display_name))
    return results
