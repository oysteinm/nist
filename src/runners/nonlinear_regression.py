"""Nonlinear Least Squares (NLS) runner for all 27 NIST StRD datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable
import warnings

import mpmath as mp
import numpy as np
import scipy.linalg
from scipy.optimize import least_squares

from src.lre import compute_lre, min_lre
from src.parser import NonlinearRegressionDataset, parse_nls_dat

PI_NIST = 3.141592653589793238462643383279

# Order matching the Stata 17 NIST StRD certification report
NLS_ORDER = [
    ("Misra1a", "Misra 1a"),
    ("Chwirut2", "Chwirut 2"),
    ("Chwirut1", "Chwirut 1"),
    ("Lanczos3", "Lanczos 3"),
    ("Gauss1", "Gauss 1"),
    ("Gauss2", "Gauss 2"),
    ("DanWood", "Daniel Wood"),
    ("Misra1b", "Misra 1b"),
    ("Kirby2", "Kirby 2"),
    ("Hahn1", "Hahn 1"),
    ("Nelson", "Nelson"),
    ("MGH17", "MGH 17"),
    ("Lanczos1", "Lanczos 1"),
    ("Lanczos2", "Lanczos 2"),
    ("Gauss3", "Gauss 3"),
    ("Misra1c", "Misra 1c"),
    ("Misra1d", "Misra 1d"),
    ("Roszman1", "Roszman 1"),
    ("ENSO", "ENSO"),
    ("MGH09", "MGH 09"),
    ("Thurber", "Thurber"),
    ("BoxBOD", "BoxBOD"),
    ("Rat42", "Ratkowsky 2"),
    ("MGH10", "MGH 10"),
    ("Eckerle4", "Eckerle4"),
    ("Rat43", "Ratkowsky 3"),
    ("Bennett5", "Bennett 5"),
]


def _c_arctan(z: np.ndarray) -> np.ndarray:
    """Complex-step-safe arctan without branch-cut discontinuities for |Re(z)| > 1."""
    if np.iscomplexobj(z):
        re = np.real(z)
        im = np.imag(z)
        return np.arctan(re) + 1j * (im / (1.0 + re * re))
    return np.arctan(z)


MODEL_FUNCTIONS: dict[str, Callable[[np.ndarray, np.ndarray], np.ndarray]] = {
    "Misra1a": lambda b, X: b[0] * (1.0 - np.exp(-b[1] * X[:, 0])),
    "Chwirut2": lambda b, X: np.exp(-b[0] * X[:, 0]) / (b[1] + b[2] * X[:, 0]),
    "Chwirut1": lambda b, X: np.exp(-b[0] * X[:, 0]) / (b[1] + b[2] * X[:, 0]),
    "Lanczos3": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-b[3] * X[:, 0])
        + b[4] * np.exp(-b[5] * X[:, 0])
    ),
    "Gauss1": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-((X[:, 0] - b[3]) ** 2) / (b[4] ** 2))
        + b[5] * np.exp(-((X[:, 0] - b[6]) ** 2) / (b[7] ** 2))
    ),
    "Gauss2": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-((X[:, 0] - b[3]) ** 2) / (b[4] ** 2))
        + b[5] * np.exp(-((X[:, 0] - b[6]) ** 2) / (b[7] ** 2))
    ),
    "DanWood": lambda b, X: b[0] * (X[:, 0] ** b[1]),
    "Misra1b": lambda b, X: b[0] * (1.0 - (1.0 + 0.5 * b[1] * X[:, 0]) ** (-2.0)),
    "Kirby2": lambda b, X: (
        (b[0] + b[1] * X[:, 0] + b[2] * (X[:, 0] ** 2))
        / (1.0 + b[3] * X[:, 0] + b[4] * (X[:, 0] ** 2))
    ),
    "Hahn1": lambda b, X: (
        (b[0] + b[1] * X[:, 0] + b[2] * (X[:, 0] ** 2) + b[3] * (X[:, 0] ** 3))
        / (1.0 + b[4] * X[:, 0] + b[5] * (X[:, 0] ** 2) + b[6] * (X[:, 0] ** 3))
    ),
    "Nelson": lambda b, X: b[0] - b[1] * X[:, 0] * np.exp(-b[2] * X[:, 1]),
    "MGH17": lambda b, X: b[0] + b[1] * np.exp(-X[:, 0] * b[3]) + b[2] * np.exp(-X[:, 0] * b[4]),
    "Lanczos1": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-b[3] * X[:, 0])
        + b[4] * np.exp(-b[5] * X[:, 0])
    ),
    "Lanczos2": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-b[3] * X[:, 0])
        + b[4] * np.exp(-b[5] * X[:, 0])
    ),
    "Gauss3": lambda b, X: (
        b[0] * np.exp(-b[1] * X[:, 0])
        + b[2] * np.exp(-((X[:, 0] - b[3]) ** 2) / (b[4] ** 2))
        + b[5] * np.exp(-((X[:, 0] - b[6]) ** 2) / (b[7] ** 2))
    ),
    "Misra1c": lambda b, X: b[0] * (1.0 - (1.0 + 2.0 * b[1] * X[:, 0]) ** (-0.5)),
    "Misra1d": lambda b, X: b[0] * b[1] * X[:, 0] * ((1.0 + b[1] * X[:, 0]) ** (-1.0)),
    "Roszman1": lambda b, X: (
        b[0] - b[1] * X[:, 0] - _c_arctan(b[2] / (X[:, 0] - b[3])) / PI_NIST
    ),
    "ENSO": lambda b, X: (
        b[0]
        + b[1] * np.cos(2.0 * PI_NIST * X[:, 0] / 12.0)
        + b[2] * np.sin(2.0 * PI_NIST * X[:, 0] / 12.0)
        + b[4] * np.cos(2.0 * PI_NIST * X[:, 0] / b[3])
        + b[5] * np.sin(2.0 * PI_NIST * X[:, 0] / b[3])
        + b[7] * np.cos(2.0 * PI_NIST * X[:, 0] / b[6])
        + b[8] * np.sin(2.0 * PI_NIST * X[:, 0] / b[6])
    ),
    "MGH09": lambda b, X: (
        b[0] * (X[:, 0] ** 2 + X[:, 0] * b[1]) / (X[:, 0] ** 2 + X[:, 0] * b[2] + b[3])
    ),
    "Thurber": lambda b, X: (
        (b[0] + b[1] * X[:, 0] + b[2] * (X[:, 0] ** 2) + b[3] * (X[:, 0] ** 3))
        / (1.0 + b[4] * X[:, 0] + b[5] * (X[:, 0] ** 2) + b[6] * (X[:, 0] ** 3))
    ),
    "BoxBOD": lambda b, X: b[0] * (1.0 - np.exp(-b[1] * X[:, 0])),
    "Rat42": lambda b, X: b[0] / (1.0 + np.exp(b[1] - b[2] * X[:, 0])),
    "MGH10": lambda b, X: b[0] * np.exp(b[1] / (X[:, 0] + b[2])),
    "Eckerle4": lambda b, X: (b[0] / b[1]) * np.exp(-0.5 * (((X[:, 0] - b[2]) / b[1]) ** 2)),
    "Rat43": lambda b, X: b[0] / ((1.0 + np.exp(b[1] - b[2] * X[:, 0])) ** (1.0 / b[3])),
    "Bennett5": lambda b, X: b[0] * ((b[1] + X[:, 0]) ** (-1.0 / b[2])),
}


def _eval_mp_point(name: str, b: list[Any], x_row: list[Any], pi_mp: Any) -> Any:
    """Evaluate the scalar model prediction in mpmath arbitrary precision for one observation."""
    x = x_row[0]
    if name in {"Misra1a", "BoxBOD"}:
        return b[0] * (1 - mp.exp(-b[1] * x))
    if name in {"Chwirut1", "Chwirut2"}:
        return mp.exp(-b[0] * x) / (b[1] + b[2] * x)
    if name in {"Lanczos1", "Lanczos2", "Lanczos3"}:
        return b[0] * mp.exp(-b[1] * x) + b[2] * mp.exp(-b[3] * x) + b[4] * mp.exp(-b[5] * x)
    if name in {"Gauss1", "Gauss2", "Gauss3"}:
        return (
            b[0] * mp.exp(-b[1] * x)
            + b[2] * mp.exp(-((x - b[3]) ** 2) / (b[4] ** 2))
            + b[5] * mp.exp(-((x - b[6]) ** 2) / (b[7] ** 2))
        )
    if name == "DanWood":
        return b[0] * (x ** b[1])
    if name == "Misra1b":
        return b[0] * (1 - (1 + b[1] * x / 2) ** (-2))
    if name == "Misra1c":
        return b[0] * (1 - (1 + 2 * b[1] * x) ** (-mp.mpf("0.5")))
    if name == "Misra1d":
        return b[0] * b[1] * x / (1 + b[1] * x)
    if name == "Kirby2":
        return (b[0] + b[1] * x + b[2] * (x**2)) / (1 + b[3] * x + b[4] * (x**2))
    if name in {"Hahn1", "Thurber"}:
        return (b[0] + b[1] * x + b[2] * (x**2) + b[3] * (x**3)) / (
            1 + b[4] * x + b[5] * (x**2) + b[6] * (x**3)
        )
    if name == "Nelson":
        return b[0] - b[1] * x_row[0] * mp.exp(-b[2] * x_row[1])
    if name == "MGH17":
        return b[0] + b[1] * mp.exp(-x * b[3]) + b[2] * mp.exp(-x * b[4])
    if name == "Roszman1":
        return b[0] - b[1] * x - mp.atan(b[2] / (x - b[3])) / pi_mp
    if name == "ENSO":
        return (
            b[0]
            + b[1] * mp.cos(2 * pi_mp * x / 12)
            + b[2] * mp.sin(2 * pi_mp * x / 12)
            + b[4] * mp.cos(2 * pi_mp * x / b[3])
            + b[5] * mp.sin(2 * pi_mp * x / b[3])
            + b[7] * mp.cos(2 * pi_mp * x / b[6])
            + b[8] * mp.sin(2 * pi_mp * x / b[6])
        )
    if name == "MGH09":
        return b[0] * (x**2 + x * b[1]) / (x**2 + x * b[2] + b[3])
    if name == "Rat42":
        return b[0] / (1 + mp.exp(b[1] - b[2] * x))
    if name == "MGH10":
        return b[0] * mp.exp(b[1] / (x + b[2]))
    if name == "Eckerle4":
        return (b[0] / b[1]) * mp.exp(-mp.mpf("0.5") * (((x - b[2]) / b[1]) ** 2))
    if name == "Rat43":
        return b[0] / ((1 + mp.exp(b[1] - b[2] * x)) ** (1 / b[3]))
    if name == "Bennett5":
        return b[0] * ((b[1] + x) ** (-1 / b[2]))
    raise KeyError(f"Unknown NLS model {name}")


def _run_mpmath_nls(ds: NonlinearRegressionDataset, init_coef: list[float], dps: int = 50) -> dict[str, float]:
    """Perform 50-digit arbitrary-precision Gauss-Newton refinement in mpmath from raw .dat strings."""
    with mp.workdps(dps):
        pi_mp = mp.mpf("3.141592653589793238462643383279")
        n = len(ds.raw_rows)
        p = len(init_coef)
        df_r = n - p

        y_mp = [mp.log(mp.mpf(r[0])) if ds.name == "Nelson" else mp.mpf(r[0]) for r in ds.raw_rows]
        X_mp = [[mp.mpf(tok) for tok in r[1:]] for r in ds.raw_rows]
        b = [mp.mpf(str(v)) for v in init_coef]
        h = mp.mpf("1e-20")

        for _ in range(5):
            r_vec = mp.matrix(n, 1)
            J_mat = mp.matrix(n, p)
            for i in range(n):
                f0 = _eval_mp_point(ds.name, b, X_mp[i], pi_mp)
                r_vec[i, 0] = y_mp[i] - f0
                for j in range(p):
                    b_plus = list(b)
                    b_plus[j] += h
                    b_minus = list(b)
                    b_minus[j] -= h
                    df = (_eval_mp_point(ds.name, b_plus, X_mp[i], pi_mp) - _eval_mp_point(ds.name, b_minus, X_mp[i], pi_mp)) / (2 * h)
                    J_mat[i, j] = -df

            jtj = J_mat.T * J_mat
            jtr = -J_mat.T * r_vec
            delta = mp.lu_solve(jtj, jtr)
            for j in range(p):
                b[j] += delta[j, 0]
            if max(abs(delta[j, 0]) for j in range(p)) < mp.mpf("1e-30"):
                break

        # Final residuals and Jacobian at 50-digit solution
        r_vec = mp.matrix(n, 1)
        J_mat = mp.matrix(n, p)
        for i in range(n):
            f0 = _eval_mp_point(ds.name, b, X_mp[i], pi_mp)
            r_vec[i, 0] = y_mp[i] - f0
            for j in range(p):
                b_plus = list(b)
                b_plus[j] += h
                b_minus = list(b)
                b_minus[j] -= h
                df = (_eval_mp_point(ds.name, b_plus, X_mp[i], pi_mp) - _eval_mp_point(ds.name, b_minus, X_mp[i], pi_mp)) / (2 * h)
                J_mat[i, j] = -df

        rss_mp = sum(r_vec[i, 0] ** 2 for i in range(n))
        jtj_inv = (J_mat.T * J_mat) ** -1
        se_mp = [mp.sqrt(jtj_inv[j, j] * (rss_mp / df_r)) for j in range(p)]

        calc_coef = [float(b[j]) for j in range(p)]
        calc_se = [float(se_mp[j]) for j in range(p)]
        calc_rss = float(rss_mp)

        return {
            "lre_coef": min_lre(calc_coef, ds.certified_coef) or 0.0,
            "lre_se": min_lre(calc_se, ds.certified_se) or 0.0,
            "lre_rss": compute_lre(calc_rss, ds.certified_rss) or 0.0,
        }


@dataclass
class NlsFitAttempt:
    start_label: str
    success: bool
    calc_coef: list[float]
    calc_se: list[float]
    calc_rss: float
    lre_coef: float
    lre_se: float
    lre_rss: float
    param_lre_coef: list[float]
    param_lre_se: list[float]


@dataclass
class NonlinearRegressionResult:
    dataset_name: str
    display_name: str
    difficulty: str
    n_obs: int
    n_params: int
    param_names: list[str]
    used_start2: bool
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
    start1_attempt: NlsFitAttempt | None = None
    start2_attempt: NlsFitAttempt | None = None
    fd_comparison: dict[str, float | bool] = field(default_factory=dict)
    alt_methods: dict[str, dict[str, float | bool]] = field(default_factory=dict)


def _complex_step_jacobian(
    f_model: Callable[[np.ndarray, np.ndarray], np.ndarray],
    b: np.ndarray,
    X: np.ndarray,
    h: float = 1e-20,
) -> np.ndarray:
    """Compute the exact Jacobian of residuals r(b) = y - f(b, X) via complex-step differentiation."""
    n = X.shape[0]
    p = len(b)
    J = np.empty((n, p), dtype=np.float64)
    for j in range(p):
        b_c = b.astype(np.complex128)
        b_c[j] += 1j * h
        f_c = f_model(b_c, X)
        # Residual is y - f(b, X), so dr/db_j = -df/db_j
        J[:, j] = -np.imag(f_c) / h
    return J


def _compute_se_from_jacobian(J: np.ndarray, rss: float, df_r: int) -> np.ndarray:
    """Compute parameter standard errors from the residual Jacobian J using pivoted QR."""
    p = J.shape[1]
    try:
        _, R, P = scipy.linalg.qr(J, mode="economic", pivoting=True)
        R_inv = scipy.linalg.solve_triangular(R, np.eye(p))
        var_piv = np.sum(R_inv**2, axis=1) * (rss / df_r)
        se = np.empty(p, dtype=np.float64)
        se[P] = np.sqrt(var_piv)
        return se
    except Exception:
        try:
            cov = scipy.linalg.pinv(J.T @ J) * (rss / df_r)
            return np.sqrt(np.maximum(0.0, np.diag(cov)))
        except Exception:
            return np.full(p, np.nan)


def _fit_from_starting_values(
    ds: NonlinearRegressionDataset,
    f_model: Callable[[np.ndarray, np.ndarray], np.ndarray],
    y_target: np.ndarray,
    start_vals: list[float],
    start_label: str,
    use_exact_jac: bool = True,
    method: str = "lm",
    jac_mode: str = "2-point",
) -> NlsFitAttempt:
    b0 = np.array(start_vals, dtype=np.float64)
    X = ds.X
    df_r = len(y_target) - len(b0)

    def residuals(b: np.ndarray) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred = f_model(b, X)
            res = y_target - pred
            if not np.all(np.isfinite(res)):
                return np.full_like(y_target, 1e100)
            return res

    def jac(b: np.ndarray) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            J = _complex_step_jacobian(f_model, b, X)
            if not np.all(np.isfinite(J)):
                return np.zeros((len(y_target), len(b)), dtype=np.float64)
            return J

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opt = least_squares(
                residuals,
                b0,
                jac=jac if use_exact_jac else jac_mode,
                method=method,
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
                max_nfev=20000,
            )
        calc_coef_arr = opt.x
        final_res = y_target - f_model(calc_coef_arr, X)
        calc_rss = float(np.dot(final_res, final_res))
        J_final = _complex_step_jacobian(f_model, calc_coef_arr, X) if use_exact_jac else opt.jac
        calc_se_arr = _compute_se_from_jacobian(J_final, calc_rss, df_r)

        calc_coef = [float(v) for v in calc_coef_arr]
        calc_se = [float(v) for v in calc_se_arr]
        param_lre_coef = [compute_lre(c, t) or 0.0 for c, t in zip(calc_coef, ds.certified_coef)]
        param_lre_se = [compute_lre(c, t) or 0.0 for c, t in zip(calc_se, ds.certified_se)]
        lre_coef = min_lre(calc_coef, ds.certified_coef) or 0.0
        lre_se = min_lre(calc_se, ds.certified_se) or 0.0
        lre_rss = compute_lre(calc_rss, ds.certified_rss) or 0.0
        converged = bool(opt.status > 0) and lre_coef >= 2.0
    except Exception:
        calc_coef = [float("nan")] * len(b0)
        calc_se = [float("nan")] * len(b0)
        calc_rss = float("nan")
        param_lre_coef = [0.0] * len(b0)
        param_lre_se = [0.0] * len(b0)
        lre_coef = 0.0
        lre_se = 0.0
        lre_rss = 0.0
        converged = False

    return NlsFitAttempt(
        start_label=start_label,
        success=converged,
        calc_coef=calc_coef,
        calc_se=calc_se,
        calc_rss=calc_rss,
        lre_coef=lre_coef,
        lre_se=lre_se,
        lre_rss=lre_rss,
        param_lre_coef=param_lre_coef,
        param_lre_se=param_lre_se,
    )


def run_nls_dataset(ds: NonlinearRegressionDataset, display_name: str | None = None) -> NonlinearRegressionResult:
    """Run NLS on a single NIST dataset, preserving the primary run and recording alternative solvers."""
    f_model = MODEL_FUNCTIONS[ds.name]
    y_target = np.log(ds.y) if ds.name == "Nelson" else ds.y

    # Primary solver (Initial run preserved): LM + Complex-Step Exact Jacobian
    attempt1 = _fit_from_starting_values(ds, f_model, y_target, ds.start1, "Start 1", use_exact_jac=True)
    attempt2 = _fit_from_starting_values(ds, f_model, y_target, ds.start2, "Start 2", use_exact_jac=True)

    if attempt1.success:
        chosen = attempt1
        used_start2 = False
    else:
        chosen = attempt2
        used_start2 = True

    # Alternative 1: LM with 2-point finite-difference Jacobian (attempt Start 1, fallback to Start 2)
    fd1 = _fit_from_starting_values(ds, f_model, y_target, ds.start1, "FD Start 1", use_exact_jac=False, method="lm", jac_mode="2-point")
    fd2 = _fit_from_starting_values(ds, f_model, y_target, ds.start2, "FD Start 2", use_exact_jac=False, method="lm", jac_mode="2-point")
    fd_chosen = fd1 if fd1.success else fd2
    fd_used_s2 = not fd1.success

    # Alternative 2: TRF with 3-point central-difference Jacobian (attempt Start 1, fallback to Start 2)
    trf1 = _fit_from_starting_values(ds, f_model, y_target, ds.start1, "TRF Start 1", use_exact_jac=False, method="trf", jac_mode="3-point")
    trf2 = _fit_from_starting_values(ds, f_model, y_target, ds.start2, "TRF Start 2", use_exact_jac=False, method="trf", jac_mode="3-point")
    trf_chosen = trf1 if trf1.success else trf2
    trf_used_s2 = not trf1.success

    # Alternative 3: Extended Precision (50-digit mpmath Gauss-Newton refinement)
    mp_res = _run_mpmath_nls(ds, chosen.calc_coef, dps=50)

    fd_summary: dict[str, float | bool] = {
        "lre_coef": fd_chosen.lre_coef,
        "lre_se": fd_chosen.lre_se,
        "lre_rss": fd_chosen.lre_rss,
        "used_start2": fd_used_s2,
    }

    alt_methods: dict[str, dict[str, float | bool]] = {
        "lm_2pt_fd": fd_summary,
        "trf_3pt_fd": {
            "lre_coef": trf_chosen.lre_coef,
            "lre_se": trf_chosen.lre_se,
            "lre_rss": trf_chosen.lre_rss,
            "used_start2": trf_used_s2,
        },
        "mpmath_50digit": {
            "lre_coef": mp_res["lre_coef"],
            "lre_se": mp_res["lre_se"],
            "lre_rss": mp_res["lre_rss"],
            "used_start2": False,
        },
    }

    return NonlinearRegressionResult(
        dataset_name=ds.name,
        display_name=display_name or ds.name,
        difficulty=ds.difficulty,
        n_obs=len(y_target),
        n_params=len(ds.param_names),
        param_names=ds.param_names,
        used_start2=used_start2,
        calc_coef=chosen.calc_coef,
        cert_coef=ds.certified_coef,
        param_lre_coef=chosen.param_lre_coef,
        lre_coef=chosen.lre_coef,
        calc_se=chosen.calc_se,
        cert_se=ds.certified_se,
        param_lre_se=chosen.param_lre_se,
        lre_se=chosen.lre_se,
        calc_rss=chosen.calc_rss,
        cert_rss=ds.certified_rss,
        lre_rss=chosen.lre_rss,
        start1_attempt=attempt1,
        start2_attempt=attempt2,
        fd_comparison=fd_summary,
        alt_methods=alt_methods,
    )


def run_nls_suite(nist_nls_dir: Path) -> list[NonlinearRegressionResult]:
    """Run all 27 NIST StRD Nonlinear Least Squares datasets in Stata report order."""
    results: list[NonlinearRegressionResult] = []
    for file_stem, display_name in NLS_ORDER:
        dat_path = nist_nls_dir / f"{file_stem}.dat"
        ds = parse_nls_dat(dat_path)
        results.append(run_nls_dataset(ds, display_name=display_name))
    return results
