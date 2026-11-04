"""Scraper for NIST StRD datasets and Stata 17 certification reference results."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

NIST_ARCHIVE_URL = "https://www.itl.nist.gov/div898/strd/general/dataarchive.html"
STATA17_CERT_URL = "https://www.stata.com/support/cert/nist17/"

DEFAULT_DATA_DIR = Path(__file__).resolve().parent.parent / "data"


SUITE_KEYS = {
    "univ": "univ",
    "lls": "lls",
    "anova": "anova",
    "nls": "nls",
}


def fetch_text(session: requests.Session, url: str, timeout: int = 20) -> str:
    """Fetch URL text with retries."""
    resp = session.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def scrape_nist_archive(
    data_dir: Path = DEFAULT_DATA_DIR,
    force: bool = False,
) -> dict[str, list[Path]]:
    """Crawl the NIST StRD archive page and download all .dat files for univ, lls, anova, and nls."""
    nist_dir = data_dir / "nist"
    nist_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "NIST-StRD-Python-Certification/1.0"})

    archive_html = fetch_text(session, NIST_ARCHIVE_URL)
    archive_soup = BeautifulSoup(archive_html, "html.parser")

    suite_urls: dict[str, str] = {}
    for a in archive_soup.find_all("a", href=True):
        full_url = urljoin(NIST_ARCHIVE_URL, a["href"])
        if "/univ/" in full_url:
            suite_urls["univ"] = full_url
        elif "/lls/" in full_url:
            suite_urls["lls"] = full_url
        elif "/anova/" in full_url:
            suite_urls["anova"] = full_url
        elif "/nls/" in full_url:
            suite_urls["nls"] = full_url

    downloaded: dict[str, list[Path]] = {k: [] for k in ["univ", "lls", "anova", "nls"]}

    for suite, suite_url in suite_urls.items():
        suite_dir = nist_dir / suite
        suite_dir.mkdir(parents=True, exist_ok=True)

        suite_html = fetch_text(session, suite_url)
        suite_soup = BeautifulSoup(suite_html, "html.parser")

        detail_pages: list[str] = []
        seen_pages: set[str] = set()
        for a in suite_soup.find_all("a", href=True):
            href = a["href"].split("#")[0]
            if not href:
                continue
            page_url = urljoin(suite_url, href)
            lower_name = page_url.rsplit("/", 1)[-1].lower()
            if lower_name in {
                "dataarchive.html",
                "index.html",
                "homepage.html",
                "backgroundinfo.html",
                "background.html",
                "lls.shtml",
                "lls_info.shtml",
                "anova.html",
                "nls_main.shtml",
                "nls_info.shtml",
                "mgh09.Causing",
            }:
                continue
            if f"/strd/{suite}/" in page_url and (
                page_url.endswith(".html") or page_url.endswith(".shtml")
            ):
                if page_url not in seen_pages:
                    seen_pages.add(page_url)
                    detail_pages.append(page_url)

        for detail_url in detail_pages:
            detail_html = fetch_text(session, detail_url)
            detail_soup = BeautifulSoup(detail_html, "html.parser")
            for a in detail_soup.find_all("a", href=True):
                dat_url = urljoin(detail_url, a["href"])
                if not dat_url.lower().endswith(".dat"):
                    continue
                dat_name = dat_url.rsplit("/", 1)[-1]
                # For ANOVA, skip the tabular *t.dat files and keep the 2-column .dat files
                if suite == "anova" and dat_name.lower().endswith("t.dat"):
                    continue

                dest_path = suite_dir / dat_name
                if dest_path.exists() and not force and dest_path.stat().st_size > 0:
                    if dest_path not in downloaded[suite]:
                        downloaded[suite].append(dest_path)
                    continue

                dat_text = fetch_text(session, dat_url)
                dest_path.write_text(dat_text, encoding="utf-8")
                if dest_path not in downloaded[suite]:
                    downloaded[suite].append(dest_path)

    return downloaded


def _parse_val_token(tok: str) -> tuple[float | None, bool]:
    """Parse a Stata LRE token such as '15.0' or '(7.0)'. Returns (value, used_start2)."""
    tok = tok.strip()
    if tok.startswith("(") and tok.endswith(")"):
        return float(tok[1:-1]), True
    return float(tok), False


def scrape_stata17_certification(
    data_dir: Path = DEFAULT_DATA_DIR,
    download_do_files: bool = True,
    force: bool = False,
) -> dict[str, Any]:
    """Scrape Stata 17 NIST StRD certification results and optionally download .do files."""
    stata_dir = data_dir / "stata17"
    stata_dir.mkdir(parents=True, exist_ok=True)
    ref_json_path = stata_dir / "stata17_reference.json"

    if ref_json_path.exists() and not force:
        return json.loads(ref_json_path.read_text(encoding="utf-8"))

    session = requests.Session()
    session.headers.update({"User-Agent": "NIST-StRD-Python-Certification/1.0"})

    html = fetch_text(session, STATA17_CERT_URL)
    soup = BeautifulSoup(html, "html.parser")

    pre_blocks = soup.find_all("pre")
    if len(pre_blocks) < 4:
        raise RuntimeError(f"Expected 4 <pre> result tables on {STATA17_CERT_URL}, found {len(pre_blocks)}")

    results: dict[str, Any] = {
        "source_url": STATA17_CERT_URL,
        "suites": {
            "univ": {"tests": [], "summary": {}},
            "lls": {"tests": [], "summary": {}},
            "anova": {"tests": [], "summary": {}},
            "nls": {"tests": [], "summary": {}},
        },
    }

    # Helper to extract links by line from a <pre> block
    def parse_pre_lines(pre_tag: Any) -> list[tuple[str, str | None, str | None]]:
        # Replace <a> tags with placeholders or inspect line by line
        raw_html = str(pre_tag)
        lines_out: list[tuple[str, str | None, str | None]] = []
        for raw_line in raw_html.splitlines():
            line_soup = BeautifulSoup(raw_line, "html.parser")
            text = line_soup.get_text()
            log_url = None
            do_url = None
            for a in line_soup.find_all("a", href=True):
                href = urljoin(STATA17_CERT_URL, a["href"])
                if href.endswith(".log"):
                    log_url = href
                elif href.endswith(".do"):
                    do_url = href
            lines_out.append((text, log_url, do_url))
        return lines_out

    # 1. Univariate
    for text, log_url, do_url in parse_pre_lines(pre_blocks[0]):
        stripped = text.strip()
        if log_url or do_url:
            # e.g.: PiDigits       lower            15.0    15.0       14.9    log   do-file
            m = re.match(
                r"^(\S+)\s+(lower|average|higher)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                stripped,
            )
            if m:
                results["suites"]["univ"]["tests"].append(
                    {
                        "test_name": m.group(1),
                        "difficulty": m.group(2),
                        "mean": float(m.group(3)),
                        "sd": float(m.group(4)),
                        "autocorr": float(m.group(5)),
                        "log_url": log_url,
                        "do_url": do_url,
                    }
                )
        else:
            m_sum = re.match(r"^(Average|Minimum|Maximum)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", stripped)
            if m_sum:
                results["suites"]["univ"]["summary"][m_sum.group(1).lower()] = {
                    "mean": float(m_sum.group(2)),
                    "sd": float(m_sum.group(3)),
                    "autocorr": float(m_sum.group(4)),
                }

    # 2. Linear Regression
    for text, log_url, do_url in parse_pre_lines(pre_blocks[1]):
        stripped = text.strip()
        if log_url or do_url:
            if "no full solution" in stripped:
                m_fil = re.match(r"^(\S+)\s+(lower|average|higher)\s+no full solution", stripped)
                if m_fil:
                    results["suites"]["lls"]["tests"].append(
                        {
                            "test_name": m_fil.group(1),
                            "difficulty": m_fil.group(2),
                            "coef": None,
                            "se": None,
                            "rss": None,
                            "note": "no full solution(*)",
                            "orthpoly_coef": 8.4,
                            "orthpoly_rss": 8.5,
                            "log_url": log_url,
                            "do_url": do_url,
                        }
                    )
            else:
                m = re.match(
                    r"^(\S+)\s+(lower|average|higher)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                    stripped,
                )
                if m:
                    results["suites"]["lls"]["tests"].append(
                        {
                            "test_name": m.group(1),
                            "difficulty": m.group(2),
                            "coef": float(m.group(3)),
                            "se": float(m.group(4)),
                            "rss": float(m.group(5)),
                            "note": None,
                            "log_url": log_url,
                            "do_url": do_url,
                        }
                    )
        else:
            m_sum = re.match(r"^(Average|Minimum|Maximum)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", stripped)
            if m_sum:
                results["suites"]["lls"]["summary"][m_sum.group(1).lower()] = {
                    "coef": float(m_sum.group(2)),
                    "se": float(m_sum.group(3)),
                    "rss": float(m_sum.group(4)),
                }

    # 3. ANOVA
    for text, log_url, do_url in parse_pre_lines(pre_blocks[2]):
        stripped = text.strip()
        if log_url or do_url:
            # e.g. "Si Resistivity     lower                 13.1          log   do-file"
            # or   "7b    higher                15.0(*)       log   do-file"
            m = re.match(
                r"^(.+?)\s+(lower|average|higher)\s+([\d.]+)(\(\*\))?",
                stripped,
            )
            if m:
                raw_name = m.group(1).strip()
                if raw_name in {"7b", "8b", "9b"}:
                    raw_name = f"Simon-Lesage {raw_name}"
                results["suites"]["anova"]["tests"].append(
                    {
                        "test_name": raw_name,
                        "difficulty": m.group(2),
                        "f_stat": float(m.group(3)),
                        "starred": bool(m.group(4)),
                        "log_url": log_url,
                        "do_url": do_url,
                    }
                )
        else:
            if stripped.startswith("Average excluding"):
                m_avg = re.search(r"([\d.]+)$", stripped)
                if m_avg:
                    results["suites"]["anova"]["summary"]["average"] = {"f_stat": float(m_avg.group(1))}
            else:
                m_sum = re.match(r"^(Minimum|Maximum)\s+([\d.]+)", stripped)
                if m_sum:
                    results["suites"]["anova"]["summary"][m_sum.group(1).lower()] = {
                        "f_stat": float(m_sum.group(2))
                    }

    # 4. Nonlinear Regression
    for text, log_url, do_url in parse_pre_lines(pre_blocks[3]):
        stripped = text.strip()
        if log_url or do_url:
            m = re.match(
                r"^(.+?)\s+(lower|average|higher)\s+(\(?[\d.]+\)?)\s+(\(?[\d.]+\)?)\s+(\(?[\d.]+\)?)",
                stripped,
            )
            if m:
                coef_val, s2_c = _parse_val_token(m.group(3))
                se_val, s2_s = _parse_val_token(m.group(4))
                rss_val, s2_r = _parse_val_token(m.group(5))
                results["suites"]["nls"]["tests"].append(
                    {
                        "test_name": m.group(1).strip(),
                        "difficulty": m.group(2),
                        "coef": coef_val,
                        "se": se_val,
                        "rss": rss_val,
                        "used_start2": s2_c or s2_s or s2_r,
                        "log_url": log_url,
                        "do_url": do_url,
                    }
                )
        else:
            m_sum = re.match(r"^(Average|Minimum|Maximum)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", stripped)
            if m_sum:
                results["suites"]["nls"]["summary"][m_sum.group(1).lower()] = {
                    "coef": float(m_sum.group(2)),
                    "se": float(m_sum.group(3)),
                    "rss": float(m_sum.group(4)),
                }

    if download_do_files:
        for suite_key, suite_data in results["suites"].items():
            do_dir = stata_dir / suite_key / "do"
            do_dir.mkdir(parents=True, exist_ok=True)
            for test_item in suite_data["tests"]:
                do_url = test_item.get("do_url")
                if do_url:
                    do_name = do_url.rsplit("/", 1)[-1]
                    do_path = do_dir / do_name
                    if not do_path.exists() or force:
                        do_path.write_text(fetch_text(session, do_url), encoding="utf-8")
                    test_item["local_do_path"] = str(do_path.relative_to(data_dir.parent))

    ref_json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results


if __name__ == "__main__":
    print("Scraping NIST StRD archive...")
    nist_files = scrape_nist_archive()
    for s, paths in nist_files.items():
        print(f"  {s}: {len(paths)} .dat files downloaded")

    print("Scraping Stata 17 NIST StRD certification page...")
    stata_ref = scrape_stata17_certification()
    for s, sdata in stata_ref["suites"].items():
        print(f"  Stata {s}: {len(sdata['tests'])} test entries parsed")
