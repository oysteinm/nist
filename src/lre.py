"""Log Relative Error (LRE) calculation utilities following Stata's NIST StRD methodology."""

from __future__ import annotations

import math
from typing import Iterable


def compute_lre(calculated: float | None, certified: float | None, max_digits: float = 15.0) -> float | None:
    """Compute Decimal Digits of Accuracy via Log Relative Error (LRE).

    Formula (per Stata 17 NIST StRD certification report):
      - LRE = min(15, -log10(|c - t| / |t|))  if |t| != 0
      - LRE = min(15, -log10(|c - t|))        otherwise
    Floored at 0.0 for valid finite numbers. Returns None if calculated is None/NaN/Inf.
    """
    if calculated is None or certified is None:
        return None
    if not (math.isfinite(calculated) and math.isfinite(certified)):
        return None

    diff = abs(calculated - certified)
    if diff == 0.0:
        return max_digits

    if certified != 0.0:
        err = diff / abs(certified)
    else:
        err = diff

    if err <= 0.0:
        return max_digits

    lre = min(max_digits, -math.log10(err))
    return max(0.0, lre)


def min_lre(calculated_values: Iterable[float | None], certified_values: Iterable[float]) -> float | None:
    """Compute the minimum LRE across multiple parameters (e.g. regression coefficients or S.E.s).

    If any parameter is None/NaN (e.g., dropped due to collinearity), returns None.
    """
    calc_list = list(calculated_values)
    cert_list = list(certified_values)
    if not calc_list or len(calc_list) != len(cert_list):
        return None

    lres: list[float] = []
    for c, t in zip(calc_list, cert_list):
        val = compute_lre(c, t)
        if val is None:
            return None
        lres.append(val)
    return min(lres) if lres else None
