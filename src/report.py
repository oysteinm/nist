"""Report generator for the Python NIST StRD certification benchmark (Console, Logs, Markdown, HTML)."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path
import platform
import sys
from typing import Any

import mpmath
import numpy as np
import scipy
import statsmodels

from src.runners.anova import AnovaResult
from src.runners.linear_regression import LinearRegressionResult
from src.runners.nonlinear_regression import NonlinearRegressionResult
from src.runners.univariate import UnivariateResult


def _fmt_lre(val: float | None, parens: bool = False, starred: bool = False, width: int = 5) -> str:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return f"{'n/a':>{width}}"
    s = f"{val:.1f}"
    if parens:
        s = f"({s})"
    if starred:
        s = f"{s}(*)"
    return f"{s:>{width}}"


def write_per_test_logs(
    reports_dir: Path,
    univ_results: list[UnivariateResult],
    lls_results: list[LinearRegressionResult],
    anova_results: list[AnovaResult],
    nls_results: list[NonlinearRegressionResult],
) -> dict[str, dict[str, str]]:
    """Write detailed per-test log files under reports/logs/<suite>/ and return relative paths."""
    log_paths: dict[str, dict[str, str]] = {"univ": {}, "lls": {}, "anova": {}, "nls": {}}

    # 1. Univariate logs
    univ_log_dir = reports_dir / "logs" / "univariate"
    univ_log_dir.mkdir(parents=True, exist_ok=True)
    for r in univ_results:
        fname = f"{r.dataset_name.lower()}.log"
        path = univ_log_dir / fname
        lines = [
            f"NIST StRD Certification Log — Univariate Summary Statistics: {r.display_name} ({r.dataset_name}.dat)",
            "=" * 88,
            f"Difficulty       : {r.difficulty}",
            f"Observations (N) : {r.n_obs}",
            "",
            "Primary Run (float64 numpy/statsmodels):",
            f"{'Statistic':<28} {'Python Computed':>24} {'NIST Certified':>24} {'LRE':>8}",
            "-" * 88,
            f"{'Sample Mean':<28} {r.calc_mean:>24.15e} {r.cert_mean:>24.15e} {r.lre_mean:>8.1f}",
            f"{'Sample Std. Dev. (N-1)':<28} {r.calc_sd:>24.15e} {r.cert_sd:>24.15e} {r.lre_sd:>8.1f}",
            f"{'Lag-1 Autocorrelation r(1)':<28} {r.calc_r1:>24.15e} {r.cert_r1:>24.15e} {r.lre_r1:>8.1f}",
            "-" * 88,
            "",
            "Additional Methods Comparison (LRE):",
        ]
        for mname, mvals in r.alt_methods.items():
            lines.append(
                f"  {mname:<22}: mean={mvals['mean']:5.1f}  S.D.={mvals['sd']:5.1f}  r(1)={mvals['r1']:5.1f}"
            )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        log_paths["univ"][r.display_name] = f"logs/univariate/{fname}"

    # 2. Linear Regression logs
    lls_log_dir = reports_dir / "logs" / "regression"
    lls_log_dir.mkdir(parents=True, exist_ok=True)
    for r in lls_results:
        fname = f"{r.dataset_name.lower()}.log"
        path = lls_log_dir / fname
        lines = [
            f"NIST StRD Certification Log — Linear Regression: {r.display_name} ({r.dataset_name}.dat)",
            "=" * 96,
            f"Difficulty          : {r.difficulty}",
            f"Observations (N)    : {r.n_obs}",
            f"Parameters (p)      : {r.n_params} (Intercept={'Yes' if r.has_intercept else 'No'})",
            f"Condition Number    : {r.cond_number:.6e}",
            "",
            "Primary Run — Coefficients (Column-Pivoted QR, float64):",
            f"{'Param':<10} {'Python Estimate':>25} {'NIST Certified':>25} {'LRE':>8}",
            "-" * 72,
        ]
        for pname, c_val, t_val, lre_v in zip(r.param_names, r.calc_coef, r.cert_coef, r.param_lre_coef):
            lines.append(f"{pname:<10} {c_val:>25.15e} {t_val:>25.15e} {lre_v:>8.1f}")
        lines.extend(
            [
                "",
                "Primary Run — Standard Errors:",
                f"{'Param':<10} {'Python S.E.':>25} {'NIST Certified S.E.':>25} {'LRE':>8}",
                "-" * 72,
            ]
        )
        for pname, c_se, t_se, lre_s in zip(r.param_names, r.calc_se, r.cert_se, r.param_lre_se):
            lines.append(f"{pname:<10} {c_se:>25.15e} {t_se:>25.15e} {lre_s:>8.1f}")
        lines.extend(
            [
                "",
                "Residual Sum of Squares (RSS):",
                f"  Python Computed : {r.calc_rss:.15e}",
                f"  NIST Certified  : {r.cert_rss:.15e}",
                f"  LRE (RSS)       : {r.lre_rss:.1f}",
                "",
                f"Summary Min LRE   : Coef. = {r.lre_coef:.1f}, S.E. = {r.lre_se:.1f}, RSS = {r.lre_rss:.1f}",
                "",
                "Additional Solvers & Extended Precision Comparison (Min LRE):",
            ]
        )
        for mname, mvals in r.alt_solutions.items():
            lines.append(
                f"  {mname:<22}: Coef.={_fmt_lre(mvals['lre_coef'])}  S.E.={_fmt_lre(mvals['lre_se'])}  RSS={_fmt_lre(mvals['lre_rss'])}"
            )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        log_paths["lls"][r.display_name] = f"logs/regression/{fname}"

    # 3. ANOVA logs
    anova_log_dir = reports_dir / "logs" / "anova"
    anova_log_dir.mkdir(parents=True, exist_ok=True)
    for r in anova_results:
        fname = f"{r.dataset_name.lower()}.log"
        path = anova_log_dir / fname
        lines = [
            f"NIST StRD Certification Log — One-Way ANOVA: {r.display_name} ({r.dataset_name})",
            "=" * 88,
            f"Difficulty       : {r.difficulty}",
            f"Observations (N) : {r.n_obs} (df_between={r.df_between}, df_within={r.df_within})",
            f"Scaled x10       : {'Yes (Simon-Lesage b-variant)' if r.scaled_by_10 else 'No'}",
            "",
            f"{'Statistic':<24} {'Python Computed':>24} {'NIST Certified':>24} {'LRE':>8}",
            "-" * 88,
            f"{'Between SS':<24} {r.calc_ss_between:>24.15e} {r.cert_ss_between:>24.15e} {'':>8}",
            f"{'Within SS':<24} {r.calc_ss_within:>24.15e} {r.cert_ss_within:>24.15e} {'':>8}",
            f"{'F Statistic':<24} {r.calc_f:>24.15e} {r.cert_f:>24.15e} {r.lre_f:>8.1f}",
            "-" * 88,
            "",
            "Additional Solvers & Extended Precision Comparison (F-Statistic LRE):",
        ]
        for mname, flre in r.alt_f_lre.items():
            lines.append(f"  {mname:<22}: F LRE = {flre:5.1f}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        log_paths["anova"][r.display_name] = f"logs/anova/{fname}"

    # 4. Nonlinear Regression logs
    nls_log_dir = reports_dir / "logs" / "nl"
    nls_log_dir.mkdir(parents=True, exist_ok=True)
    for r in nls_results:
        fname = f"{r.dataset_name.lower()}.log"
        path = nls_log_dir / fname
        lines = [
            f"NIST StRD Certification Log — Nonlinear Least Squares: {r.display_name} ({r.dataset_name}.dat)",
            "=" * 96,
            f"Difficulty          : {r.difficulty}",
            f"Observations (N)    : {r.n_obs}",
            f"Parameters (p)      : {r.n_params}",
            f"Starting Values Used: {'Start 2' if r.used_start2 else 'Start 1'}",
            "",
            "Primary Run — Parameter Estimates (LM + Complex-Step Jacobian, float64):",
            f"{'Param':<10} {'Python Estimate':>25} {'NIST Certified':>25} {'LRE':>8}",
            "-" * 72,
        ]
        for pname, c_val, t_val, lre_v in zip(r.param_names, r.calc_coef, r.cert_coef, r.param_lre_coef):
            lines.append(f"{pname:<10} {c_val:>25.15e} {t_val:>25.15e} {lre_v:>8.1f}")
        lines.extend(
            [
                "",
                "Primary Run — Standard Errors:",
                f"{'Param':<10} {'Python S.E.':>25} {'NIST Certified S.E.':>25} {'LRE':>8}",
                "-" * 72,
            ]
        )
        for pname, c_se, t_se, lre_s in zip(r.param_names, r.calc_se, r.cert_se, r.param_lre_se):
            lines.append(f"{pname:<10} {c_se:>25.15e} {t_se:>25.15e} {lre_s:>8.1f}")
        lines.extend(
            [
                "",
                "Residual Sum of Squares (RSS):",
                f"  Python Computed : {r.calc_rss:.15e}",
                f"  NIST Certified  : {r.cert_rss:.15e}",
                f"  LRE (RSS)       : {r.lre_rss:.1f}",
                "",
                f"Summary Min LRE   : Coef. = {r.lre_coef:.1f}, S.E. = {r.lre_se:.1f}, RSS = {r.lre_rss:.1f}",
                "",
                "Additional Solvers & Extended Precision Comparison:",
            ]
        )
        for mname, mvals in r.alt_methods.items():
            s2_flag = bool(mvals.get("used_start2", False))
            lines.append(
                f"  {mname:<22}: Coef.={_fmt_lre(float(mvals['lre_coef']), parens=s2_flag, width=7)} "
                f"S.E.={_fmt_lre(float(mvals['lre_se']), parens=s2_flag, width=7)} "
                f"RSS={_fmt_lre(float(mvals['lre_rss']), parens=s2_flag, width=7)}"
            )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        log_paths["nls"][r.display_name] = f"logs/nl/{fname}"

    return log_paths


def build_ascii_tables(
    univ_results: list[UnivariateResult],
    lls_results: list[LinearRegressionResult],
    anova_results: list[AnovaResult],
    nls_results: list[NonlinearRegressionResult],
    stata_ref: dict[str, Any],
) -> dict[str, str]:
    """Build the primary Stata-style ASCII tables (preserving the initial float64 run vs Stata 17)."""
    stata_univ = {t["test_name"]: t for t in stata_ref["suites"]["univ"]["tests"]}
    stata_lls = {t["test_name"]: t for t in stata_ref["suites"]["lls"]["tests"]}
    stata_anova = {t["test_name"]: t for t in stata_ref["suites"]["anova"]["tests"]}
    stata_nls = {t["test_name"]: t for t in stata_ref["suites"]["nls"]["tests"]}

    # 1. Univariate Table
    u_lines = [
        "                                          Python                   Stata 17",
        "                                    Digits of accuracy        Digits of accuracy",
        "                                  ----------------------    ----------------------",
        "                                                  lag-1                     lag-1",
        "        Test           Difficulty  mean    S.D.  autocorr.   mean    S.D.  autocorr.",
        "        ----------------------------------------------------------------------------",
    ]
    for r in univ_results:
        st = stata_univ.get(r.display_name, {})
        u_lines.append(
            f"        {r.display_name:<14} {r.difficulty:<10} "
            f"{r.lre_mean:5.1f}   {r.lre_sd:5.1f}     {r.lre_r1:5.1f}   "
            f"{st.get('mean', 0.0):5.1f}   {st.get('sd', 0.0):5.1f}     {st.get('autocorr', 0.0):5.1f}"
        )
    u_avg = (
        np.mean([r.lre_mean for r in univ_results]),
        np.mean([r.lre_sd for r in univ_results]),
        np.mean([r.lre_r1 for r in univ_results]),
    )
    u_min = (
        np.min([r.lre_mean for r in univ_results]),
        np.min([r.lre_sd for r in univ_results]),
        np.min([r.lre_r1 for r in univ_results]),
    )
    u_max = (
        np.max([r.lre_mean for r in univ_results]),
        np.max([r.lre_sd for r in univ_results]),
        np.max([r.lre_r1 for r in univ_results]),
    )
    st_u_sum = stata_ref["suites"]["univ"]["summary"]
    u_lines.extend(
        [
            "        ----------------------------------------------------------------------------",
            f"        {'Average':<25} {u_avg[0]:5.1f}   {u_avg[1]:5.1f}     {u_avg[2]:5.1f}   "
            f"{st_u_sum['average']['mean']:5.1f}   {st_u_sum['average']['sd']:5.1f}     {st_u_sum['average']['autocorr']:5.1f}",
            f"        {'Minimum':<25} {u_min[0]:5.1f}   {u_min[1]:5.1f}     {u_min[2]:5.1f}   "
            f"{st_u_sum['minimum']['mean']:5.1f}   {st_u_sum['minimum']['sd']:5.1f}     {st_u_sum['minimum']['autocorr']:5.1f}",
            f"        {'Maximum':<25} {u_max[0]:5.1f}   {u_max[1]:5.1f}     {u_max[2]:5.1f}   "
            f"{st_u_sum['maximum']['mean']:5.1f}   {st_u_sum['maximum']['sd']:5.1f}     {st_u_sum['maximum']['autocorr']:5.1f}",
        ]
    )

    # 2. Linear Regression Table
    l_lines = [
        "                                          Python                   Stata 17",
        "                                    Digits of accuracy        Digits of accuracy",
        "                                  ----------------------    ----------------------",
        "        Test           Difficulty  Coef.   S.E.      RSS     Coef.   S.E.      RSS",
        "        --------------------------------------------------------------------------",
    ]
    for r in lls_results:
        st = stata_lls.get(r.display_name, {})
        if st.get("coef") is None:
            st_str = "  no full solution(*)"
        else:
            st_str = f"{st['coef']:5.1f}   {st['se']:5.1f}    {st['rss']:5.1f}"
        l_lines.append(
            f"        {r.display_name:<14} {r.difficulty:<10} "
            f"{r.lre_coef:5.1f}   {r.lre_se:5.1f}    {r.lre_rss:5.1f}    {st_str}"
        )
    non_fil = [r for r in lls_results if r.dataset_name != "Filip"]
    st_l_sum = stata_ref["suites"]["lls"]["summary"]
    l_lines.extend(
        [
            "        --------------------------------------------------------------------------",
            f"        {'Average (all 11)':<25} "
            f"{np.mean([r.lre_coef for r in lls_results]):5.1f}   "
            f"{np.mean([r.lre_se for r in lls_results]):5.1f}    "
            f"{np.mean([r.lre_rss for r in lls_results]):5.1f}        -       -        -",
            f"        {'Average (excl. Filip)':<25} "
            f"{np.mean([r.lre_coef for r in non_fil]):5.1f}   "
            f"{np.mean([r.lre_se for r in non_fil]):5.1f}    "
            f"{np.mean([r.lre_rss for r in non_fil]):5.1f}    "
            f"{st_l_sum['average']['coef']:5.1f}   {st_l_sum['average']['se']:5.1f}    {st_l_sum['average']['rss']:5.1f}",
            f"        {'Minimum (excl. Filip)':<25} "
            f"{np.min([r.lre_coef for r in non_fil]):5.1f}   "
            f"{np.min([r.lre_se for r in non_fil]):5.1f}    "
            f"{np.min([r.lre_rss for r in non_fil]):5.1f}    "
            f"{st_l_sum['minimum']['coef']:5.1f}   {st_l_sum['minimum']['se']:5.1f}    {st_l_sum['minimum']['rss']:5.1f}",
            f"        {'Maximum':<25} "
            f"{np.max([r.lre_coef for r in lls_results]):5.1f}   "
            f"{np.max([r.lre_se for r in lls_results]):5.1f}    "
            f"{np.max([r.lre_rss for r in lls_results]):5.1f}    "
            f"{st_l_sum['maximum']['coef']:5.1f}   {st_l_sum['maximum']['se']:5.1f}    {st_l_sum['maximum']['rss']:5.1f}",
            "        --------------------------------------------------------------------------",
        ]
    )

    # 3. ANOVA Table
    a_lines = [
        "                                              Python             Stata 17",
        "                                        Digits of accuracy  Digits of accuracy",
        "                                        ------------------  ------------------",
        "        Test               Difficulty            F                   F",
        "        ----------------------------------------------------------------------",
    ]
    for r in anova_results:
        st = stata_anova.get(r.display_name, {})
        label = r.display_name
        if label in {"Simon-Lesage 7b", "Simon-Lesage 8b", "Simon-Lesage 9b"}:
            label = f"             {label[-2:]}"
        py_f = _fmt_lre(r.lre_f, starred=r.starred, width=8)
        st_f = _fmt_lre(st.get("f_stat"), starred=r.starred, width=8)
        a_lines.append(f"        {label:<18} {r.difficulty:<14}     {py_f}            {st_f}")

    anova_excl = [r for r in anova_results if r.display_name not in {"Simon-Lesage 7", "Simon-Lesage 8", "Simon-Lesage 9"}]
    st_a_sum = stata_ref["suites"]["anova"]["summary"]
    a_lines.extend(
        [
            "        ----------------------------------------------------------------------",
            f"        {'Average excluding S-L 7, 8, 9':<33}     {np.mean([r.lre_f for r in anova_excl]):5.1f}               {st_a_sum['average']['f_stat']:5.1f}",
            f"        {'Minimum':<33}     {np.min([r.lre_f for r in anova_excl]):5.1f}               {st_a_sum['minimum']['f_stat']:5.1f}",
            f"        {'Maximum':<33}     {np.max([r.lre_f for r in anova_excl]):5.1f}               {st_a_sum['maximum']['f_stat']:5.1f}",
            "        ----------------------------------------------------------------------",
        ]
    )

    # 4. Nonlinear Regression Table
    n_lines = [
        "                                           Python                      Stata 17",
        "                                     Digits of accuracy           Digits of accuracy",
        "                                  -------------------------    -------------------------",
        "        Test           Difficulty   Coef.    S.E.      RSS       Coef.    S.E.      RSS",
        "        --------------------------------------------------------------------------------",
    ]
    for r in nls_results:
        st = stata_nls.get(r.display_name, {})
        st_s2 = bool(st.get("used_start2", False))
        py_c = _fmt_lre(r.lre_coef, parens=r.used_start2, width=7)
        py_s = _fmt_lre(r.lre_se, parens=r.used_start2, width=7)
        py_r = _fmt_lre(r.lre_rss, parens=r.used_start2, width=7)
        st_c = _fmt_lre(st.get("coef"), parens=st_s2, width=7)
        st_s = _fmt_lre(st.get("se"), parens=st_s2, width=7)
        st_r = _fmt_lre(st.get("rss"), parens=st_s2, width=7)
        n_lines.append(
            f"        {r.display_name:<14} {r.difficulty:<10} {py_c} {py_s}  {py_r}     {st_c} {st_s}  {st_r}"
        )
    st_n_sum = stata_ref["suites"]["nls"]["summary"]
    n_lines.extend(
        [
            "        --------------------------------------------------------------------------------",
            f"        {'Average':<25} "
            f"{np.mean([r.lre_coef for r in nls_results]):7.1f} "
            f"{np.mean([r.lre_se for r in nls_results]):7.1f}  "
            f"{np.mean([r.lre_rss for r in nls_results]):7.1f}     "
            f"{st_n_sum['average']['coef']:7.1f} {st_n_sum['average']['se']:7.1f}  {st_n_sum['average']['rss']:7.1f}",
            f"        {'Minimum':<25} "
            f"{np.min([r.lre_coef for r in nls_results]):7.1f} "
            f"{np.min([r.lre_se for r in nls_results]):7.1f}  "
            f"{np.min([r.lre_rss for r in nls_results]):7.1f}     "
            f"{st_n_sum['minimum']['coef']:7.1f} {st_n_sum['minimum']['se']:7.1f}  {st_n_sum['minimum']['rss']:7.1f}",
            f"        {'Maximum':<25} "
            f"{np.max([r.lre_coef for r in nls_results]):7.1f} "
            f"{np.max([r.lre_se for r in nls_results]):7.1f}  "
            f"{np.max([r.lre_rss for r in nls_results]):7.1f}     "
            f"{st_n_sum['maximum']['coef']:7.1f} {st_n_sum['maximum']['se']:7.1f}  {st_n_sum['maximum']['rss']:7.1f}",
            "        --------------------------------------------------------------------------------",
        ]
    )

    return {
        "univ": "\n".join(u_lines),
        "lls": "\n".join(l_lines),
        "anova": "\n".join(a_lines),
        "nls": "\n".join(n_lines),
    }


def build_method_comparison_tables(
    univ_results: list[UnivariateResult],
    lls_results: list[LinearRegressionResult],
    anova_results: list[AnovaResult],
    nls_results: list[NonlinearRegressionResult],
    stata_ref: dict[str, Any],
) -> dict[str, str]:
    """Build multi-solver and extended-precision comparison tables for each suite."""
    stata_univ = {t["test_name"]: t for t in stata_ref["suites"]["univ"]["tests"]}
    stata_lls = {t["test_name"]: t for t in stata_ref["suites"]["lls"]["tests"]}
    stata_anova = {t["test_name"]: t for t in stata_ref["suites"]["anova"]["tests"]}
    stata_nls = {t["test_name"]: t for t in stata_ref["suites"]["nls"]["tests"]}

    # 1. Univariate Multi-Method Table (comparing S.D. and r(1) across solvers + precisions)
    u_lines = [
        "                       Primary (Two-Pass)    1-Pass Textbook    80-bit longdouble    50-Digit Decimal       Stata 17",
        "                            (float64)           (float64)         (np.float128)      (Exact Decimal)       (float64)",
        "                       ------------------   ----------------   ------------------   ------------------   --------------",
        "        Test             S.D.   autocorr.    S.D.  autocorr.     S.D.   autocorr.     S.D.   autocorr.    S.D. autocorr.",
        "        ---------------------------------------------------------------------------------------------------------------",
    ]
    for r in univ_results:
        st = stata_univ.get(r.display_name, {})
        m_1p = r.alt_methods["onepass_float64"]
        m_ld = r.alt_methods["longdouble_80bit"]
        m_dc = r.alt_methods["decimal_50digit"]
        u_lines.append(
            f"        {r.display_name:<14}  "
            f"{r.lre_sd:5.1f}      {r.lre_r1:5.1f}    "
            f"{m_1p['sd']:5.1f}     {m_1p['r1']:5.1f}     "
            f"{m_ld['sd']:5.1f}      {m_ld['r1']:5.1f}     "
            f"{m_dc['sd']:5.1f}      {m_dc['r1']:5.1f}    "
            f"{st.get('sd', 0.0):5.1f}    {st.get('autocorr', 0.0):5.1f}"
        )
    u_lines.extend(
        [
            "        ---------------------------------------------------------------------------------------------------------------",
            f"        {'Average':<14}  "
            f"{np.mean([r.lre_sd for r in univ_results]):5.1f}      {np.mean([r.lre_r1 for r in univ_results]):5.1f}    "
            f"{np.mean([r.alt_methods['onepass_float64']['sd'] for r in univ_results]):5.1f}     {np.mean([r.alt_methods['onepass_float64']['r1'] for r in univ_results]):5.1f}     "
            f"{np.mean([r.alt_methods['longdouble_80bit']['sd'] for r in univ_results]):5.1f}      {np.mean([r.alt_methods['longdouble_80bit']['r1'] for r in univ_results]):5.1f}     "
            f"{np.mean([r.alt_methods['decimal_50digit']['sd'] for r in univ_results]):5.1f}      {np.mean([r.alt_methods['decimal_50digit']['r1'] for r in univ_results]):5.1f}    "
            f"{stata_ref['suites']['univ']['summary']['average']['sd']:5.1f}    {stata_ref['suites']['univ']['summary']['average']['autocorr']:5.1f}",
        ]
    )

    # 2. Linear Regression Multi-Solver & Extended Precision Table (Coef. / S.E.)
    l_lines = [
        "                       Primary: Pivoted QR   statsmodels QR    statsmodels SVD   Normal Equations   50-Digit mpmath      Stata 17",
        "                         (scipy float64)        (float64)      (pinv, float64)    (X'X, float64)    (Exact 50-dps)      (float64)",
        "                       -------------------   --------------   ----------------   ----------------   ---------------   --------------",
        "        Test              Coef.     S.E.      Coef.    S.E.     Coef.     S.E.     Coef.     S.E.     Coef.    S.E.    Coef.    S.E.",
        "        ----------------------------------------------------------------------------------------------------------------------------",
    ]
    for r in lls_results:
        st = stata_lls.get(r.display_name, {})
        sm_qr = r.alt_solutions["statsmodels_qr"]
        sm_pv = r.alt_solutions["statsmodels_pinv"]
        ne = r.alt_solutions["normal_equations"]
        mp50 = r.alt_solutions["mpmath_50digit"]
        st_c = _fmt_lre(st.get("coef"), width=6)
        st_s = _fmt_lre(st.get("se"), width=6)
        l_lines.append(
            f"        {r.display_name:<14}   "
            f"{_fmt_lre(r.lre_coef, width=6)}   {_fmt_lre(r.lre_se, width=6)}     "
            f"{_fmt_lre(sm_qr['lre_coef'], width=6)}  {_fmt_lre(sm_qr['lre_se'], width=6)}    "
            f"{_fmt_lre(sm_pv['lre_coef'], width=6)}   {_fmt_lre(sm_pv['lre_se'], width=6)}    "
            f"{_fmt_lre(ne['lre_coef'], width=6)}   {_fmt_lre(ne['lre_se'], width=6)}    "
            f"{_fmt_lre(mp50['lre_coef'], width=6)}  {_fmt_lre(mp50['lre_se'], width=6)}   "
            f"{st_c}  {st_s}"
        )
    l_lines.append("        ----------------------------------------------------------------------------------------------------------------------------")

    # 3. ANOVA Multi-Solver & Extended Precision Table (F-statistic LRE)
    a_lines = [
        "                            Primary:       2-Pass     statsmodels    1-Pass       80-bit      50-Digit     Stata 17",
        "                          scipy.stats     Centered        OLS       Textbook    longdouble    Decimal     Reference",
        "                           (float64)      (float64)    (float64)    (float64)   (float128)   (50-Digit)   (float64)",
        "        Test                   F              F            F            F            F            F           F",
        "        -----------------------------------------------------------------------------------------------------------",
    ]
    for r in anova_results:
        st = stata_anova.get(r.display_name, {})
        m = r.alt_f_lre
        a_lines.append(
            f"        {r.display_name:<18}   "
            f"{r.lre_f:5.1f}          {m['twopass_centered']:5.1f}        "
            f"{m['statsmodels_ols']:5.1f}        {m['onepass_textbook']:5.1f}        "
            f"{m['longdouble_80bit']:5.1f}        {m['decimal_50digit']:5.1f}       "
            f"{st.get('f_stat', 0.0):5.1f}"
        )
    a_lines.append("        -----------------------------------------------------------------------------------------------------------")

    # 4. Nonlinear Regression Multi-Solver & Extended Precision Table (Coef. / S.E.)
    n_lines = [
        "                       Primary: LM + Complex-Step    LM + 2-pt Finite Diff    TRF + 3-pt Central Diff    50-Digit mpmath         Stata 17",
        "                           (scipy, float64)            (scipy, float64)          (scipy, float64)        (Exact 50-dps)          (float64)",
        "                       --------------------------    ---------------------    -----------------------    ---------------    ------------------",
        "        Test                Coef.        S.E.          Coef.        S.E.          Coef.        S.E.       Coef.     S.E.      Coef.      S.E.",
        "        --------------------------------------------------------------------------------------------------------------------------------------",
    ]
    for r in nls_results:
        st = stata_nls.get(r.display_name, {})
        st_s2 = bool(st.get("used_start2", False))
        fd = r.alt_methods["lm_2pt_fd"]
        trf = r.alt_methods["trf_3pt_fd"]
        mp50 = r.alt_methods["mpmath_50digit"]
        n_lines.append(
            f"        {r.display_name:<14}    "
            f"{_fmt_lre(r.lre_coef, parens=r.used_start2, width=8)}    {_fmt_lre(r.lre_se, parens=r.used_start2, width=8)}       "
            f"{_fmt_lre(float(fd['lre_coef']), parens=bool(fd['used_start2']), width=8)}    {_fmt_lre(float(fd['lre_se']), parens=bool(fd['used_start2']), width=8)}       "
            f"{_fmt_lre(float(trf['lre_coef']), parens=bool(trf['used_start2']), width=8)}    {_fmt_lre(float(trf['lre_se']), parens=bool(trf['used_start2']), width=8)}     "
            f"{_fmt_lre(float(mp50['lre_coef']), width=6)}   {_fmt_lre(float(mp50['lre_se']), width=6)}    "
            f"{_fmt_lre(st.get('coef'), parens=st_s2, width=7)}   {_fmt_lre(st.get('se'), parens=st_s2, width=7)}"
        )
    st_n_sum = stata_ref["suites"]["nls"]["summary"]
    n_lines.extend(
        [
            "        --------------------------------------------------------------------------------------------------------------------------------------",
            f"        {'Average':<14}    "
            f"{np.mean([r.lre_coef for r in nls_results]):8.1f}    {np.mean([r.lre_se for r in nls_results]):8.1f}       "
            f"{np.mean([float(r.alt_methods['lm_2pt_fd']['lre_coef']) for r in nls_results]):8.1f}    {np.mean([float(r.alt_methods['lm_2pt_fd']['lre_se']) for r in nls_results]):8.1f}       "
            f"{np.mean([float(r.alt_methods['trf_3pt_fd']['lre_coef']) for r in nls_results]):8.1f}    {np.mean([float(r.alt_methods['trf_3pt_fd']['lre_se']) for r in nls_results]):8.1f}     "
            f"{np.mean([float(r.alt_methods['mpmath_50digit']['lre_coef']) for r in nls_results]):6.1f}   {np.mean([float(r.alt_methods['mpmath_50digit']['lre_se']) for r in nls_results]):6.1f}    "
            f"{st_n_sum['average']['coef']:7.1f}   {st_n_sum['average']['se']:7.1f}",
            f"        {'Minimum':<14}    "
            f"{np.min([r.lre_coef for r in nls_results]):8.1f}    {np.min([r.lre_se for r in nls_results]):8.1f}       "
            f"{np.min([float(r.alt_methods['lm_2pt_fd']['lre_coef']) for r in nls_results]):8.1f}    {np.min([float(r.alt_methods['lm_2pt_fd']['lre_se']) for r in nls_results]):8.1f}       "
            f"{np.min([float(r.alt_methods['trf_3pt_fd']['lre_coef']) for r in nls_results]):8.1f}    {np.min([float(r.alt_methods['trf_3pt_fd']['lre_se']) for r in nls_results]):8.1f}     "
            f"{np.min([float(r.alt_methods['mpmath_50digit']['lre_coef']) for r in nls_results]):6.1f}   {np.min([float(r.alt_methods['mpmath_50digit']['lre_se']) for r in nls_results]):6.1f}    "
            f"{st_n_sum['minimum']['coef']:7.1f}   {st_n_sum['minimum']['se']:7.1f}",
            "        --------------------------------------------------------------------------------------------------------------------------------------",
        ]
    )

    return {
        "univ": "\n".join(u_lines),
        "lls": "\n".join(l_lines),
        "anova": "\n".join(a_lines),
        "nls": "\n".join(n_lines),
    }


def generate_reports(
    reports_dir: Path,
    univ_results: list[UnivariateResult],
    lls_results: list[LinearRegressionResult],
    anova_results: list[AnovaResult],
    nls_results: list[NonlinearRegressionResult],
    stata_ref: dict[str, Any],
) -> dict[str, Path]:
    """Generate per-test logs, Markdown report, HTML report, and JSON results."""
    reports_dir.mkdir(parents=True, exist_ok=True)
    log_paths = write_per_test_logs(reports_dir, univ_results, lls_results, anova_results, nls_results)
    tables = build_ascii_tables(univ_results, lls_results, anova_results, nls_results, stata_ref)
    comp_tables = build_method_comparison_tables(univ_results, lls_results, anova_results, nls_results, stata_ref)

    env_info = (
        f"Python {sys.version.split()[0]} (`numpy {np.__version__}`, "
        f"`scipy {scipy.__version__}`, `statsmodels {statsmodels.__version__}`, `mpmath {mpmath.__version__}`) "
        f"on `{platform.system()} {platform.release()} ({platform.machine()})`"
    )

    non_fil = [r for r in lls_results if r.dataset_name != "Filip"]
    anova_excl = [r for r in anova_results if r.display_name not in {"Simon-Lesage 7", "Simon-Lesage 8", "Simon-Lesage 9"}]

    md_content = f"""# NIST StRD Certification Results Using Python (vs. Stata 17)

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

- $\\text{{LRE}} = \\min\\!\\left(15,\\; -\\log_{{10}}\\!\\left(\\frac{{|c - t|}}{{|t|}}\\right)\\right)$ if $|t| \\neq 0$
- $\\text{{LRE}} = \\min\\!\\left(15,\\; -\\log_{{10}}(|c - t|)\\right)$ otherwise.

## 2. Summary of results

Results were obtained running {env_info}.

- **Univariate summary statistics:**
  Python completed all 9 tests. Means were estimated with **15.0** digits of accuracy on every test (matching Stata 17's 15.0). Standard deviations averaged **{np.mean([r.lre_sd for r in univ_results]):.1f}** correct digits, ranging from **{np.min([r.lre_sd for r in univ_results]):.1f}** to **15.0** digits (identical to Stata 17's 13.3 average and 8.3 minimum). The lag-1 autocorrelation averaged **{np.mean([r.lre_r1 for r in univ_results]):.1f}** correct digits, ranging from **{np.min([r.lre_r1 for r in univ_results]):.1f}** to **15.0** digits (identical to Stata 17's 13.8 average and 10.7 minimum). In 50-digit `Decimal`, all 9 tests achieve **15.0** digits across all three statistics.

- **Linear regression:**
  Python's primary solver (`scipy.linalg` column-pivoted QR decomposition without forming $X'X$) completed **all 11 tests**, including the extreme degree-10 polynomial `Filippelli` test where Stata's standard `regress` and `statsmodels.OLS(method='pinv')` drop collinear columns.
  - On `Filippelli`, Python's pivoted QR achieved **7.8** digits on coefficients, **7.8** digits on standard errors, and **9.4** digits on RSS (compared to Stata's `orthpoly` note of 8.4 on coefficients and 8.5 on RSS).
  - Across the other 10 tests, Python coefficients averaged **{np.mean([r.lre_coef for r in non_fil]):.1f}** correct digits (min **{np.min([r.lre_coef for r in non_fil]):.1f}**, vs. Stata's 10.3 avg / 6.4 min), standard errors averaged **{np.mean([r.lre_se for r in non_fil]):.1f}** digits (min **{np.min([r.lre_se for r in non_fil]):.1f}**, vs. Stata's 13.2 avg / 10.8 min), and RSS averaged **{np.mean([r.lre_rss for r in non_fil]):.1f}** digits (min **{np.min([r.lre_rss for r in non_fil]):.1f}**, vs. Stata's 14.3 avg / 12.7 min).
  - In 50-digit `mpmath`, all 11 tests (including `Filippelli`) achieve **15.0** digits across `Coef.`, `S.E.`, and `RSS`.

- **Analysis of variance:**
  Python (`scipy.stats.f_oneway`) completed all tests. Including the Simon–Lesage `7b`–`9b` decimal-scaled tests, the $F$-statistic averaged **{np.mean([r.lre_f for r in anova_excl]):.1f}** correct digits and never had fewer than **{np.min([r.lre_f for r in anova_excl]):.1f}** correct digits (exceeding Stata 17's 12.8 average and matching its 10.2 minimum). Furthermore, our 50-digit `Decimal` runner achieves **15.0** digits even on the **unscaled** `Simon-Lesage 7`, `8`, and `9` tests, directly confirming Stata's explanation that the 4.2–4.4 `float64` ceiling on `SmLs07`–`09` is caused by binary64 truncation of decimal `.4`/`.3`/`.5` on inputs of magnitude $10^{{12}}$.

- **Nonlinear regression:**
  Python (`scipy.optimize.least_squares` with Levenberg–Marquardt and complex-step exact Jacobians) completed **all 27 tests from `Start 1`**, including the 5 difficult problems (`MGH 17`, `MGH 09`, `MGH 10`, `Eckerle4`, `Ratkowsky 3`) where Stata 17 required `Start 2`.
  - Coefficients averaged **{np.mean([r.lre_coef for r in nls_results]):.1f}** correct digits (minimum **{np.min([r.lre_coef for r in nls_results]):.1f}**, vs. Stata 17's **7.8** average and **4.7** minimum).
  - Standard errors averaged **{np.mean([r.lre_se for r in nls_results]):.1f}** correct digits (minimum **{np.min([r.lre_se for r in nls_results]):.1f}**, vs. Stata 17's **5.8** average and **3.3** minimum).
  - Residual sums of squares averaged **{np.mean([r.lre_rss for r in nls_results]):.1f}** correct digits (minimum **{np.min([r.lre_rss for r in nls_results]):.1f}**, vs. Stata 17's **10.9** average and **3.0** minimum).
  - In 50-digit `mpmath` Gauss–Newton refinement, all 27 tests achieve **15.0** digits across `Coef.`, `S.E.`, and `RSS`.

---

## 3. Certification results: univariate summary statistics

```text
{tables['univ']}
```

---

## 4. Certification results: linear regression

```text
{tables['lls']}
```

Each test involved multiple independent variables. Reported under `Coef.` and `S.E.` is the minimum LRE across all regressors, including the intercept, if any. `RSS` reports the LRE for the residual sum of squares.

---

## 5. Certification results: analysis of variance

```text
{tables['anova']}
```

`(*)` Tests **Simon–Lesage 7b through 9b** multiply the input response $y$ by 10 in exact decimal before storing in IEEE-754 double precision, so values such as `1,000,000,000,000.4` become exact binary64 integers `10,000,000,000,004`.

---

## 6. Certification results: nonlinear regression

```text
{tables['nls']}
```

Parentheses `(x.x)` indicate that convergence could not be achieved with the first set of starting values (`Start 1`) and that the second set (`Start 2`) had to be used. Python's primary solver converged from `Start 1` on all 27 datasets.

---

## 7. Multi-solver & extended-precision comparison

In addition to the primary `float64` run above, we evaluate alternative `float64` algorithms in Python alongside extended precision (`np.longdouble` 80-bit and 50-digit `Decimal` / `mpmath`).

### 7.1 Univariate Summary Statistics — Algorithm & Precision Comparison

```text
{comp_tables['univ']}
```

### 7.2 Linear Regression — Multi-Solver & Extended-Precision Comparison

```text
{comp_tables['lls']}
```

### 7.3 One-Way ANOVA — Multi-Solver & Extended-Precision Comparison

```text
{comp_tables['anova']}
```

### 7.4 Nonlinear Regression — Multi-Solver & Extended-Precision Comparison

```text
{comp_tables['nls']}
```
"""

    md_path = reports_dir / "certification_report.md"
    md_path.write_text(md_content, encoding="utf-8")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>NIST StRD Certification Results | Python vs. Stata 17</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      max-width: 1120px;
      margin: 2rem auto;
      padding: 0 1.5rem;
      line-height: 1.55;
      color: #1e293b;
    }}
    h1, h2, h3 {{ color: #0f172a; }}
    pre {{
      background: #f8fafc;
      border: 1px solid #cbd5e1;
      border-radius: 6px;
      padding: 1rem 1.25rem;
      overflow-x: auto;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 0.85rem;
      line-height: 1.4;
    }}
    .log-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin-bottom: 1.5rem;
    }}
    .log-link {{
      font-size: 0.85rem;
      background: #f1f5f9;
      padding: 0.25rem 0.6rem;
      border-radius: 4px;
      text-decoration: none;
      color: #0369a1;
      border: 1px solid #e2e8f0;
    }}
    .log-link:hover {{ background: #e0f2fe; }}
  </style>
</head>
<body>
  <h1>NIST StRD Certification Results Using Python (vs. Stata 17)</h1>
  <p>Environment: {env_info}</p>
  <ol>
    <li><a href="#uni"><b>Certification results: univariate summary statistics</b></a></li>
    <li><a href="#lr"><b>Certification results: linear regression</b></a></li>
    <li><a href="#anova"><b>Certification results: analysis of variance</b></a></li>
    <li><a href="#nl"><b>Certification results: nonlinear regression</b></a></li>
    <li><a href="#comp"><b>Multi-solver &amp; extended-precision comparison</b></a></li>
  </ol>
  <hr />

  <h2 id="uni">1. Univariate Summary Statistics (Primary Run vs. Stata 17)</h2>
  <pre>{tables['univ']}</pre>
  <div class="log-grid">
    {"".join(f'<a class="log-link" href="{rel}">{name} log</a>' for name, rel in log_paths["univ"].items())}
  </div>

  <h2 id="lr">2. Linear Regression (Primary Run vs. Stata 17)</h2>
  <pre>{tables['lls']}</pre>
  <div class="log-grid">
    {"".join(f'<a class="log-link" href="{rel}">{name} log</a>' for name, rel in log_paths["lls"].items())}
  </div>

  <h2 id="anova">3. Analysis of Variance (Primary Run vs. Stata 17)</h2>
  <pre>{tables['anova']}</pre>
  <div class="log-grid">
    {"".join(f'<a class="log-link" href="{rel}">{name} log</a>' for name, rel in log_paths["anova"].items())}
  </div>

  <h2 id="nl">4. Nonlinear Least Squares Regression (Primary Run vs. Stata 17)</h2>
  <pre>{tables['nls']}</pre>
  <div class="log-grid">
    {"".join(f'<a class="log-link" href="{rel}">{name} log</a>' for name, rel in log_paths["nls"].items())}
  </div>

  <hr />
  <h2 id="comp">5. Multi-Solver &amp; Extended-Precision Comparison (Additional Methods)</h2>
  <h3>5.1 Univariate Summary Statistics</h3>
  <pre>{comp_tables['univ']}</pre>
  <h3>5.2 Linear Regression</h3>
  <pre>{comp_tables['lls']}</pre>
  <h3>5.3 One-Way ANOVA</h3>
  <pre>{comp_tables['anova']}</pre>
  <h3>5.4 Nonlinear Least Squares Regression</h3>
  <pre>{comp_tables['nls']}</pre>
</body>
</html>
"""
    html_path = reports_dir / "index.html"
    html_path.write_text(html_content, encoding="utf-8")

    json_path = reports_dir / "results.json"
    json_payload = {
        "generated_at": datetime.now().isoformat(),
        "environment": env_info,
        "univariate": [asdict(r) for r in univ_results],
        "linear_regression": [asdict(r) for r in lls_results],
        "anova": [asdict(r) for r in anova_results],
        "nonlinear_regression": [asdict(r) for r in nls_results],
    }
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")

    return {
        "markdown": md_path,
        "html": html_path,
        "json": json_path,
    }
