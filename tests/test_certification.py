"""Unit and integration tests for the NIST StRD Python certification suite."""

from __future__ import annotations

from pathlib import Path

from src.lre import compute_lre, min_lre
from src.runners.anova import run_anova_suite
from src.runners.linear_regression import run_lls_suite
from src.runners.nonlinear_regression import run_nls_suite
from src.runners.univariate import run_univariate_suite

DATA_NIST_DIR = Path(__file__).resolve().parent.parent / "data" / "nist"


def test_lre_formula() -> None:
    assert compute_lre(10.0, 10.0) == 15.0
    assert compute_lre(0.0, 0.0) == 15.0
    assert abs((compute_lre(1.001, 1.0) or 0.0) - 3.0) < 1e-12
    assert abs((compute_lre(0.001, 0.0) or 0.0) - 3.0) < 1e-12
    assert min_lre([1.0, 1.001], [1.0, 1.0]) == compute_lre(1.001, 1.0)


def test_univariate_suite() -> None:
    results = run_univariate_suite(DATA_NIST_DIR / "univ")
    assert len(results) == 9
    for r in results:
        assert r.lre_mean == 15.0
        assert round(r.lre_sd, 1) >= 8.3
        assert round(r.lre_r1, 1) >= 10.7


def test_lls_suite() -> None:
    results = run_lls_suite(DATA_NIST_DIR / "lls")
    assert len(results) == 11
    for r in results:
        assert r.lre_coef >= 6.5
        assert r.lre_se >= 7.5
        assert r.lre_rss >= 9.0


def test_anova_suite() -> None:
    results = run_anova_suite(DATA_NIST_DIR / "anova")
    assert len(results) == 14
    for r in results:
        if r.display_name in {"Simon-Lesage 7", "Simon-Lesage 8", "Simon-Lesage 9"}:
            assert r.lre_f >= 4.0
        else:
            assert r.lre_f >= 10.1


def test_nls_suite() -> None:
    results = run_nls_suite(DATA_NIST_DIR / "nls")
    assert len(results) == 27
    for r in results:
        assert not r.used_start2
        assert r.lre_coef >= 7.0
        assert r.lre_se >= 3.0
        assert r.lre_rss >= 2.8
        assert float(r.alt_methods["mpmath_50digit"]["lre_coef"]) >= 7.9
        assert float(r.alt_methods["mpmath_50digit"]["lre_se"]) >= 8.2


def test_extended_precision_methods() -> None:
    for r in run_univariate_suite(DATA_NIST_DIR / "univ"):
        assert r.alt_methods["decimal_50digit"]["mean"] >= 14.9
        assert r.alt_methods["decimal_50digit"]["sd"] >= 14.9
        assert round(r.alt_methods["decimal_50digit"]["r1"], 1) >= 14.8

    for r in run_lls_suite(DATA_NIST_DIR / "lls"):
        assert round(r.alt_solutions["mpmath_50digit"]["lre_coef"] or 0.0, 1) >= 14.3
        assert round(r.alt_solutions["mpmath_50digit"]["lre_se"] or 0.0, 1) >= 14.5
        assert round(r.alt_solutions["mpmath_50digit"]["lre_rss"] or 0.0, 1) >= 14.5

    for r in run_anova_suite(DATA_NIST_DIR / "anova"):
        assert r.alt_f_lre["decimal_50digit"] >= 13.0
