#!/usr/bin/env python3
"""CLI entry point to scrape NIST StRD & Stata 17 data, run all certification tests, and generate reports."""

from __future__ import annotations

import argparse
from pathlib import Path

from src.report import build_ascii_tables, build_method_comparison_tables, generate_reports
from src.runners.anova import run_anova_suite
from src.runners.linear_regression import run_lls_suite
from src.runners.nonlinear_regression import run_nls_suite
from src.runners.univariate import run_univariate_suite
from src.scraper import DEFAULT_DATA_DIR, scrape_nist_archive, scrape_stata17_certification


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the NIST StRD numerical accuracy certification benchmark in Python."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory for cached NIST .dat files and Stata 17 reference data (default: ./data)",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "reports",
        help="Directory for generated reports and per-test logs (default: ./reports)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Force re-scraping of NIST StRD and Stata 17 web pages even if cached locally.",
    )
    args = parser.parse_args()

    nist_dir = args.data_dir / "nist"
    stata_json = args.data_dir / "stata17" / "stata17_reference.json"

    if args.refresh or not nist_dir.exists() or not stata_json.exists():
        print("Scraping NIST StRD archive and Stata 17 certification reference...")
        scrape_nist_archive(data_dir=args.data_dir, force=args.refresh)
    stata_ref = scrape_stata17_certification(data_dir=args.data_dir, force=args.refresh)

    print("Running NIST StRD certification suites in Python...")
    univ_results = run_univariate_suite(nist_dir / "univ")
    lls_results = run_lls_suite(nist_dir / "lls")
    anova_results = run_anova_suite(nist_dir / "anova")
    nls_results = run_nls_suite(nist_dir / "nls")

    paths = generate_reports(
        reports_dir=args.reports_dir,
        univ_results=univ_results,
        lls_results=lls_results,
        anova_results=anova_results,
        nls_results=nls_results,
        stata_ref=stata_ref,
    )

    tables = build_ascii_tables(univ_results, lls_results, anova_results, nls_results, stata_ref)
    comp_tables = build_method_comparison_tables(univ_results, lls_results, anova_results, nls_results, stata_ref)

    print("\n" + "=" * 86)
    print("1. CERTIFICATION RESULTS: UNIVARIATE SUMMARY STATISTICS (PRIMARY RUN)")
    print("=" * 86)
    print(tables["univ"])

    print("\n" + "=" * 86)
    print("2. CERTIFICATION RESULTS: LINEAR REGRESSION (PRIMARY RUN)")
    print("=" * 86)
    print(tables["lls"])

    print("\n" + "=" * 86)
    print("3. CERTIFICATION RESULTS: ANALYSIS OF VARIANCE (PRIMARY RUN)")
    print("=" * 86)
    print(tables["anova"])

    print("\n" + "=" * 86)
    print("4. CERTIFICATION RESULTS: NONLINEAR REGRESSION (PRIMARY RUN)")
    print("=" * 86)
    print(tables["nls"])

    print("\n" + "=" * 118)
    print("5. MULTI-SOLVER & EXTENDED-PRECISION COMPARISON (ADDITIONAL METHODS)")
    print("=" * 118)
    print("\n--- 5.1 Univariate Summary Statistics ---")
    print(comp_tables["univ"])
    print("\n--- 5.2 Linear Regression ---")
    print(comp_tables["lls"])
    print("\n--- 5.3 Analysis of Variance ---")
    print(comp_tables["anova"])
    print("\n--- 5.4 Nonlinear Regression ---")
    print(comp_tables["nls"])

    print("\nGenerated artifacts:")
    print(f"  - Markdown Report : {paths['markdown']}")
    print(f"  - HTML Report     : {paths['html']}")
    print(f"  - JSON Results    : {paths['json']}")
    print(f"  - Per-test Logs   : {args.reports_dir / 'logs'}")


if __name__ == "__main__":
    main()
