"""Univariate Summary Statistics runner for NIST StRD datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np
from statsmodels.tsa.stattools import acf

from src.lre import compute_lre
from src.parser import UnivariateDataset, parse_univariate_dat

# Order matching the Stata 17 NIST certification report
UNIVARIATE_ORDER = [
    ("PiDigits", "PiDigits"),
    ("Lottery", "Lottery"),
    ("Lew", "Lew"),
    ("Mavro", "Mavro"),
    ("Michelso", "Michelson"),
    ("NumAcc1", "NumAcc-1"),
    ("NumAcc2", "NumAcc-2"),
    ("NumAcc3", "NumAcc-3"),
    ("NumAcc4", "NumAcc-4"),
]


@dataclass
class UnivariateResult:
    dataset_name: str
    display_name: str
    difficulty: str
    n_obs: int
    # Initial primary run (float64 numpy/statsmodels)
    calc_mean: float
    cert_mean: float
    lre_mean: float
    calc_sd: float
    cert_sd: float
    lre_sd: float
    calc_r1: float
    cert_r1: float
    lre_r1: float
    # Additional methods (multi-solver & extended precision)
    alt_methods: dict[str, dict[str, float]] = field(default_factory=dict)


def _run_decimal_univariate(ds: UnivariateDataset) -> dict[str, float]:
    """Compute univariate statistics in 50-digit exact Decimal arithmetic from raw strings."""
    getcontext().prec = 50
    vals = [Decimal(s) for s in ds.raw_y_strings]
    n = Decimal(len(vals))
    mean_d = sum(vals) / n
    diffs = [v - mean_d for v in vals]
    ss_d = sum(d * d for d in diffs)
    sd_d = (ss_d / (n - Decimal(1))).sqrt()
    num_r1 = sum(diffs[i] * diffs[i + 1] for i in range(len(vals) - 1))
    r1_d = num_r1 / ss_d

    return {
        "mean": compute_lre(float(mean_d), ds.certified_mean) or 0.0,
        "sd": compute_lre(float(sd_d), ds.certified_sd) or 0.0,
        "r1": compute_lre(float(r1_d), ds.certified_r1) or 0.0,
    }


def _run_longdouble_univariate(ds: UnivariateDataset) -> dict[str, float]:
    """Compute univariate statistics in np.longdouble (80-bit extended precision on x86_64)."""
    y_ld = np.array([ np.longdouble(s) for s in ds.raw_y_strings ], dtype=np.longdouble)
    mean_ld = np.mean(y_ld)
    diffs = y_ld - mean_ld
    ss_ld = np.sum(diffs * diffs)
    sd_ld = np.sqrt(ss_ld / (len(y_ld) - 1))
    r1_ld = np.sum(diffs[:-1] * diffs[1:]) / ss_ld
    return {
        "mean": compute_lre(float(mean_ld), ds.certified_mean) or 0.0,
        "sd": compute_lre(float(sd_ld), ds.certified_sd) or 0.0,
        "r1": compute_lre(float(r1_ld), ds.certified_r1) or 0.0,
    }


def _run_onepass_univariate(ds: UnivariateDataset) -> dict[str, float]:
    """Compute univariate SD and r1 using naive textbook one-pass formulas in float64."""
    y = ds.y
    n = len(y)
    mean_1p = float(np.sum(y) / n)
    var_1p = float((np.sum(y * y) - n * (mean_1p**2)) / (n - 1))
    sd_1p = float(np.sqrt(var_1p)) if var_1p > 0.0 else 0.0
    # Shifted/compensated two-pass vs one-pass
    diffs = y - mean_1p
    denom = float(np.sum(y * y) - n * (mean_1p**2))
    r1_1p = float(np.sum(diffs[:-1] * diffs[1:]) / denom) if denom > 0.0 else 0.0
    return {
        "mean": compute_lre(mean_1p, ds.certified_mean) or 0.0,
        "sd": compute_lre(sd_1p, ds.certified_sd) or 0.0,
        "r1": compute_lre(r1_1p, ds.certified_r1) or 0.0,
    }


def run_univariate_dataset(ds: UnivariateDataset, display_name: str | None = None) -> UnivariateResult:
    """Compute sample mean, sample SD (ddof=1), and lag-1 sample autocorrelation and evaluate LRE."""
    y = ds.y
    calc_mean = float(np.mean(y))
    calc_sd = float(np.std(y, ddof=1))
    calc_r1 = float(acf(y, nlags=1, fft=False)[1])

    lre_mean = compute_lre(calc_mean, ds.certified_mean) or 0.0
    lre_sd = compute_lre(calc_sd, ds.certified_sd) or 0.0
    lre_r1 = compute_lre(calc_r1, ds.certified_r1) or 0.0

    alt_methods = {
        "onepass_float64": _run_onepass_univariate(ds),
        "longdouble_80bit": _run_longdouble_univariate(ds),
        "decimal_50digit": _run_decimal_univariate(ds),
    }

    return UnivariateResult(
        dataset_name=ds.name,
        display_name=display_name or ds.name,
        difficulty=ds.difficulty,
        n_obs=len(y),
        calc_mean=calc_mean,
        cert_mean=ds.certified_mean,
        lre_mean=lre_mean,
        calc_sd=calc_sd,
        cert_sd=ds.certified_sd,
        lre_sd=lre_sd,
        calc_r1=calc_r1,
        cert_r1=ds.certified_r1,
        lre_r1=lre_r1,
        alt_methods=alt_methods,
    )


def run_univariate_suite(nist_univ_dir: Path) -> list[UnivariateResult]:
    """Run all 9 NIST StRD Univariate datasets in Stata report order."""
    results: list[UnivariateResult] = []
    for file_stem, display_name in UNIVARIATE_ORDER:
        dat_path = nist_univ_dir / f"{file_stem}.dat"
        ds = parse_univariate_dat(dat_path)
        results.append(run_univariate_dataset(ds, display_name=display_name))
    return results
