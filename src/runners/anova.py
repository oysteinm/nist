"""One-Way Analysis of Variance (ANOVA) runner for NIST StRD datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, getcontext
from pathlib import Path
import warnings

import numpy as np
from scipy.stats import f_oneway
import statsmodels.api as sm

from src.lre import compute_lre
from src.parser import AnovaDataset, parse_anova_dat

# Order matching the Stata 17 NIST certification report (11 base tests + 3 scaled variants 7b, 8b, 9b)
ANOVA_ORDER = [
    ("SiRstv", "Si Resistivity", False),
    ("SmLs01", "Simon-Lesage 1", False),
    ("SmLs02", "Simon-Lesage 2", False),
    ("SmLs03", "Simon-Lesage 3", False),
    ("AtmWtAg", "Ag Atomic Wt", False),
    ("SmLs04", "Simon-Lesage 4", False),
    ("SmLs05", "Simon-Lesage 5", False),
    ("SmLs06", "Simon-Lesage 6", False),
    ("SmLs07", "Simon-Lesage 7", False),
    ("SmLs07", "Simon-Lesage 7b", True),
    ("SmLs08", "Simon-Lesage 8", False),
    ("SmLs08", "Simon-Lesage 8b", True),
    ("SmLs09", "Simon-Lesage 9", False),
    ("SmLs09", "Simon-Lesage 9b", True),
]


@dataclass
class AnovaResult:
    dataset_name: str
    display_name: str
    difficulty: str
    scaled_by_10: bool
    starred: bool
    n_obs: int
    df_between: int
    df_within: int
    calc_ss_between: float
    cert_ss_between: float
    calc_ss_within: float
    cert_ss_within: float
    calc_f: float
    cert_f: float
    lre_f: float
    alt_f_lre: dict[str, float] = field(default_factory=dict)


def _run_statsmodels_anova_f(groups: np.ndarray, y: np.ndarray) -> float:
    """Compute one-way ANOVA F-statistic via statsmodels OLS on treatment indicators."""
    unique_g = np.unique(groups)
    # Build treatment indicator matrix with intercept
    X = np.column_stack([np.ones(len(y))] + [(groups == g).astype(np.float64) for g in unique_g[1:]])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = sm.OLS(y, X).fit()
        return float(res.fvalue)


def _run_onepass_anova_f(group_arrays: list[np.ndarray], y: np.ndarray) -> float:
    """Compute one-way ANOVA F-statistic via naive textbook one-pass sums of squares in float64."""
    n = len(y)
    k = len(group_arrays)
    df_b = k - 1
    df_w = n - k
    grand_sum = float(np.sum(y))
    corr_term = (grand_sum * grand_sum) / n
    total_ss = float(np.sum(y * y)) - corr_term
    between_ss = sum((float(np.sum(g)) ** 2) / len(g) for g in group_arrays) - corr_term
    within_ss = total_ss - between_ss
    if within_ss <= 0.0 or between_ss <= 0.0:
        return 0.0
    return float((between_ss / df_b) / (within_ss / df_w))


def _run_longdouble_anova_f(ds: AnovaDataset, scale_by_10: bool) -> float:
    """Compute one-way ANOVA F-statistic in 80-bit np.longdouble from raw strings."""
    if scale_by_10:
        y_ld = np.array([np.longdouble(str(Decimal(s) * 10)) for s in ds.raw_y_strings], dtype=np.longdouble)
    else:
        y_ld = np.array([np.longdouble(s) for s in ds.raw_y_strings], dtype=np.longdouble)
    unique_g = np.unique(ds.groups)
    g_arrays = [y_ld[ds.groups == g] for g in unique_g]
    k = len(unique_g)
    n = len(y_ld)
    grand_mean = np.mean(y_ld)
    c_groups = [g - grand_mean for g in g_arrays]
    g_means = [np.mean(cg) for cg in c_groups]
    ss_b = sum(len(cg) * (gm * gm) for cg, gm in zip(c_groups, g_means))
    ss_w = sum(np.sum((cg - gm) ** 2) for cg, gm in zip(c_groups, g_means))
    return float((ss_b / (k - 1)) / (ss_w / (n - k)))


def _run_decimal_anova_f(ds: AnovaDataset, scale_by_10: bool) -> float:
    """Compute one-way ANOVA F-statistic in 50-digit exact Decimal arithmetic from raw strings."""
    getcontext().prec = 50
    factor = Decimal(10) if scale_by_10 else Decimal(1)
    y_dec = [Decimal(s) * factor for s in ds.raw_y_strings]
    n = len(y_dec)
    unique_g = sorted(set(int(g) for g in ds.groups))
    k = len(unique_g)
    grand_mean = sum(y_dec) / Decimal(n)

    ss_b = Decimal(0)
    ss_w = Decimal(0)
    for g in unique_g:
        g_vals = [y_dec[i] - grand_mean for i in range(n) if int(ds.groups[i]) == g]
        n_g = Decimal(len(g_vals))
        gm = sum(g_vals) / n_g
        ss_b += n_g * gm * gm
        ss_w += sum((v - gm) * (v - gm) for v in g_vals)

    ms_b = ss_b / Decimal(k - 1)
    ms_w = ss_w / Decimal(n - k)
    return float(ms_b / ms_w)


def run_anova_dataset(
    ds: AnovaDataset,
    display_name: str | None = None,
    scale_by_10: bool = False,
) -> AnovaResult:
    """Run One-Way ANOVA and evaluate LRE for the F-statistic across primary and alternative solvers."""
    if scale_by_10:
        y = np.array([float(Decimal(s) * 10) for s in ds.raw_y_strings], dtype=np.float64)
        cert_ss_between = ds.certified_ss_between * 100.0
        cert_ss_within = ds.certified_ss_within * 100.0
    else:
        y = ds.y
        cert_ss_between = ds.certified_ss_between
        cert_ss_within = ds.certified_ss_within

    unique_groups = np.unique(ds.groups)
    group_arrays = [y[ds.groups == g] for g in unique_groups]

    k = len(unique_groups)
    n = len(y)
    df_between = k - 1
    df_within = n - k

    # Primary solver (Initial run preserved): Scipy's f_oneway
    scipy_f = float(f_oneway(*group_arrays).statistic)

    # Explicit two-pass centered ANOVA for SS_between, SS_within, and F
    grand_mean = float(np.mean(y))
    centered_groups = [g - grand_mean for g in group_arrays]
    group_means_c = [float(np.mean(cg)) for cg in centered_groups]

    ss_between = float(sum(len(cg) * (gm**2) for cg, gm in zip(centered_groups, group_means_c)))
    ss_within = float(sum(np.sum((cg - gm) ** 2) for cg, gm in zip(centered_groups, group_means_c)))

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    calc_f = float(ms_between / ms_within)

    lre_f_scipy = compute_lre(scipy_f, ds.certified_f) or 0.0
    lre_f_twopass = compute_lre(calc_f, ds.certified_f) or 0.0
    lre_f_sm = compute_lre(_run_statsmodels_anova_f(ds.groups, y), ds.certified_f) or 0.0
    lre_f_onepass = compute_lre(_run_onepass_anova_f(group_arrays, y), ds.certified_f) or 0.0
    lre_f_ld = compute_lre(_run_longdouble_anova_f(ds, scale_by_10), ds.certified_f) or 0.0
    lre_f_dec = compute_lre(_run_decimal_anova_f(ds, scale_by_10), ds.certified_f) or 0.0

    starred = ds.name in {"SmLs07", "SmLs08", "SmLs09"}

    return AnovaResult(
        dataset_name=f"{ds.name}b" if scale_by_10 else ds.name,
        display_name=display_name or ds.name,
        difficulty=ds.difficulty,
        scaled_by_10=scale_by_10,
        starred=starred,
        n_obs=n,
        df_between=df_between,
        df_within=df_within,
        calc_ss_between=ss_between,
        cert_ss_between=cert_ss_between,
        calc_ss_within=ss_within,
        cert_ss_within=cert_ss_within,
        calc_f=scipy_f,
        cert_f=ds.certified_f,
        lre_f=lre_f_scipy,
        alt_f_lre={
            "scipy_f_oneway": lre_f_scipy,
            "twopass_centered": lre_f_twopass,
            "statsmodels_ols": lre_f_sm,
            "onepass_textbook": lre_f_onepass,
            "longdouble_80bit": lre_f_ld,
            "decimal_50digit": lre_f_dec,
        },
    )


def run_anova_suite(nist_anova_dir: Path) -> list[AnovaResult]:
    """Run all 14 ANOVA tests (11 base + 3 scaled 7b-9b variants) in Stata report order."""
    results: list[AnovaResult] = []
    for file_stem, display_name, scale_by_10 in ANOVA_ORDER:
        dat_path = nist_anova_dir / f"{file_stem}.dat"
        ds = parse_anova_dat(dat_path)
        results.append(run_anova_dataset(ds, display_name=display_name, scale_by_10=scale_by_10))
    return results
