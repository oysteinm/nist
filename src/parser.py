"""Parser for NIST StRD .dat files across Univariate, Linear Regression, ANOVA, and Nonlinear Regression."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any

import numpy as np


@dataclass
class UnivariateDataset:
    name: str
    difficulty: str
    certified_mean: float
    certified_sd: float
    certified_r1: float
    y: np.ndarray
    raw_y_strings: list[str] = field(default_factory=list)
    file_path: Path | None = None


@dataclass
class LinearRegressionDataset:
    name: str
    difficulty: str
    has_intercept: bool
    param_names: list[str]
    certified_coef: list[float]
    certified_se: list[float]
    certified_rss: float
    certified_r2: float | None
    certified_res_sd: float | None
    model_expr: str
    y: np.ndarray
    X_raw: np.ndarray
    raw_rows: list[list[str]] = field(default_factory=list)
    file_path: Path | None = None


@dataclass
class AnovaDataset:
    name: str
    difficulty: str
    certified_df_between: int
    certified_df_within: int
    certified_ss_between: float
    certified_ss_within: float
    certified_ms_between: float
    certified_ms_within: float
    certified_f: float
    certified_r2: float | None
    certified_res_sd: float | None
    groups: np.ndarray
    y: np.ndarray
    raw_y_strings: list[str] = field(default_factory=list)
    file_path: Path | None = None


@dataclass
class NonlinearRegressionDataset:
    name: str
    difficulty: str
    param_names: list[str]
    start1: list[float]
    start2: list[float]
    certified_coef: list[float]
    certified_se: list[float]
    certified_rss: float
    certified_res_sd: float | None
    df: int | None
    n_obs: int | None
    model_lines: list[str]
    y: np.ndarray
    X: np.ndarray
    raw_rows: list[list[str]] = field(default_factory=list)
    file_path: Path | None = None


def _extract_difficulty(text: str) -> str:
    m = re.search(r"(Lower|Average|Higher)\s+Level\s+of\s+Difficulty", text, re.IGNORECASE)
    if m:
        return m.group(1).lower()
    # Michelso.dat omits the difficulty line in the .dat header; it is classified as 'lower' by NIST
    return "lower"


def _extract_dataset_name(text: str, fallback: str) -> str:
    m = re.search(r"Dataset Name:\s*([A-Za-z0-9_-]+)", text)
    if m:
        return m.group(1).strip()
    return fallback


def _extract_data_line_range(text: str, total_lines: int) -> tuple[int, int]:
    """Return 1-based inclusive (start_line, end_line) for the Data section."""
    m = re.search(r"Data\s*:?\s*\(?lines\s+(\d+)\s+to\s+(\d+)\)?", text, re.IGNORECASE)
    if m:
        return int(m.group(1)), int(m.group(2))
    return 61, total_lines


def parse_univariate_dat(path: Path) -> UnivariateDataset:
    """Parse a NIST StRD Univariate .dat file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    name = _extract_dataset_name(text, path.stem)
    difficulty = _extract_difficulty(text)
    start_line, end_line = _extract_data_line_range(text, len(lines))

    m_mean = re.search(r"Sample Mean.*?:\s*([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    m_sd = re.search(r"Sample Standard Deviation.*?:\s*([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    m_r1 = re.search(r"Sample Autocorrelation.*?:\s*([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    if not (m_mean and m_sd and m_r1):
        raise ValueError(f"Failed to parse certified values in {path}")

    raw_tokens: list[str] = []
    for line in lines[start_line - 1 : end_line]:
        stripped = line.strip()
        if stripped:
            raw_tokens.extend(stripped.split())

    y = np.array([float(tok) for tok in raw_tokens], dtype=np.float64)
    return UnivariateDataset(
        name=name,
        difficulty=difficulty,
        certified_mean=float(m_mean.group(1)),
        certified_sd=float(m_sd.group(1)),
        certified_r1=float(m_r1.group(1)),
        y=y,
        raw_y_strings=raw_tokens,
        file_path=path,
    )


def parse_lls_dat(path: Path) -> LinearRegressionDataset:
    """Parse a NIST StRD Linear Least Squares .dat file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    name = _extract_dataset_name(text, path.stem)
    difficulty = _extract_difficulty(text)
    start_line, end_line = _extract_data_line_range(text, len(lines))

    param_names: list[str] = []
    certified_coef: list[float] = []
    certified_se: list[float] = []

    for line in lines[: start_line - 1]:
        m_param = re.match(
            r"^\s*(B\d+)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s*$",
            line,
        )
        if m_param:
            param_names.append(m_param.group(1))
            certified_coef.append(float(m_param.group(2)))
            certified_se.append(float(m_param.group(3)))

    has_intercept = "B0" in param_names

    m_rss = re.search(
        r"^\s*Residual\s+\d+\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)",
        text,
        re.MULTILINE,
    )
    if not m_rss:
        raise ValueError(f"Could not find Residual Sum of Squares in {path}")
    certified_rss = float(m_rss.group(1))

    m_r2 = re.search(r"R-Squared\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    certified_r2 = float(m_r2.group(1)) if m_r2 else None

    m_rsd = re.search(r"Standard Deviation\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    certified_res_sd = float(m_rsd.group(1)) if m_rsd else None

    # Extract model equation lines
    model_lines: list[str] = []
    in_model = False
    for line in lines[:35]:
        if line.strip().startswith("y =") or in_model:
            if not line.strip() and in_model:
                break
            in_model = True
            model_lines.append(line.strip())
    model_expr = " ".join(model_lines)

    rows: list[list[float]] = []
    raw_rows: list[list[str]] = []
    for line in lines[start_line - 1 : end_line]:
        stripped = line.strip()
        if stripped:
            toks = stripped.split()
            raw_rows.append(toks)
            rows.append([float(tok) for tok in toks])

    mat = np.array(rows, dtype=np.float64)
    y = mat[:, 0]
    X_raw = mat[:, 1:]

    return LinearRegressionDataset(
        name=name,
        difficulty=difficulty,
        has_intercept=has_intercept,
        param_names=param_names,
        certified_coef=certified_coef,
        certified_se=certified_se,
        certified_rss=certified_rss,
        certified_r2=certified_r2,
        certified_res_sd=certified_res_sd,
        model_expr=model_expr,
        y=y,
        X_raw=X_raw,
        raw_rows=raw_rows,
        file_path=path,
    )


def parse_anova_dat(path: Path) -> AnovaDataset:
    """Parse a NIST StRD One-Way ANOVA .dat file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    name = _extract_dataset_name(text, path.stem)
    difficulty = _extract_difficulty(text)
    start_line, end_line = _extract_data_line_range(text, len(lines))

    # Match Between and Within rows in Certified Values table
    m_between = re.search(
        r"^\s*Between\s+\S*\s*(\d+)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)",
        text,
        re.MULTILINE,
    )
    m_within = re.search(
        r"^\s*Within\s+\S*\s*(\d+)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)",
        text,
        re.MULTILINE,
    )
    if not (m_between and m_within):
        raise ValueError(f"Failed to parse ANOVA certified table in {path}")

    m_r2 = re.search(r"Certified R-Squared\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)
    m_rsd = re.search(r"Standard Deviation\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)", text)

    groups_list: list[int] = []
    y_list: list[float] = []
    raw_y_strings: list[str] = []

    for line in lines[start_line - 1 : end_line]:
        stripped = line.strip()
        if stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                groups_list.append(int(float(parts[0])))
                raw_y_strings.append(parts[1])
                y_list.append(float(parts[1]))

    return AnovaDataset(
        name=name,
        difficulty=difficulty,
        certified_df_between=int(m_between.group(1)),
        certified_ss_between=float(m_between.group(2)),
        certified_ms_between=float(m_between.group(3)),
        certified_f=float(m_between.group(4)),
        certified_df_within=int(m_within.group(1)),
        certified_ss_within=float(m_within.group(2)),
        certified_ms_within=float(m_within.group(3)),
        certified_r2=float(m_r2.group(1)) if m_r2 else None,
        certified_res_sd=float(m_rsd.group(1)) if m_rsd else None,
        groups=np.array(groups_list, dtype=np.int64),
        y=np.array(y_list, dtype=np.float64),
        raw_y_strings=raw_y_strings,
        file_path=path,
    )


def parse_nls_dat(path: Path) -> NonlinearRegressionDataset:
    """Parse a NIST StRD Nonlinear Least Squares .dat file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    name = _extract_dataset_name(text, path.stem)
    difficulty = _extract_difficulty(text)
    start_line, end_line = _extract_data_line_range(text, len(lines))

    param_names: list[str] = []
    start1: list[float] = []
    start2: list[float] = []
    certified_coef: list[float] = []
    certified_se: list[float] = []

    num_pat = r"([+-]?[\d.]+(?:[Ee][+-]?\d+)?)"
    param_line_re = re.compile(
        rf"^\s*([A-Za-z0-9_]+)\s*=\s*{num_pat}\s+{num_pat}\s+{num_pat}\s+{num_pat}\s*$"
    )
    for line in lines[: start_line - 1]:
        m = param_line_re.match(line)
        if m:
            param_names.append(m.group(1))
            start1.append(float(m.group(2)))
            start2.append(float(m.group(3)))
            certified_coef.append(float(m.group(4)))
            certified_se.append(float(m.group(5)))

    m_rss = re.search(rf"Residual Sum of Squares:\s*{num_pat}", text)
    m_rsd = re.search(rf"Residual Standard Deviation:\s*{num_pat}", text)
    m_df = re.search(r"Degrees of Freedom:\s*(\d+)", text)
    m_nobs = re.search(r"Number of Observations:\s*(\d+)", text)

    if not m_rss or not param_names:
        raise ValueError(f"Failed to parse NLS certified values in {path}")

    model_lines: list[str] = []
    in_model = False
    for line in lines[:42]:
        if "Model:" in line:
            in_model = True
            continue
        if in_model:
            if "Starting values" in line or "Certified Values" in line:
                break
            if line.strip():
                model_lines.append(line.strip())

    rows: list[list[float]] = []
    raw_rows: list[list[str]] = []
    for line in lines[start_line - 1 : end_line]:
        stripped = line.strip()
        if stripped:
            toks = stripped.split()
            raw_rows.append(toks)
            rows.append([float(tok) for tok in toks])

    mat = np.array(rows, dtype=np.float64)
    y = mat[:, 0]
    X = mat[:, 1:]

    return NonlinearRegressionDataset(
        name=name,
        difficulty=difficulty,
        param_names=param_names,
        start1=start1,
        start2=start2,
        certified_coef=certified_coef,
        certified_se=certified_se,
        certified_rss=float(m_rss.group(1)),
        certified_res_sd=float(m_rsd.group(1)) if m_rsd else None,
        df=int(m_df.group(1)) if m_df else None,
        n_obs=int(m_nobs.group(1)) if m_nobs else None,
        model_lines=model_lines,
        y=y,
        X=X,
        raw_rows=raw_rows,
        file_path=path,
    )
