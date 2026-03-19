"""
OpenAlex (via pyalex) -> export EXCEL
with execution traces printed to stdout.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from pyalex import Works, config

import os
import requests
from urllib.parse import urlparse

# ---------- LOGGING HELPERS ----------
def now() -> str:
    """Return current time as a readable string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    """Print a timestamped log message."""
    print(f"{now()} - {message}")
# ------------------------------------


# --------- SETTINGS (edit these) ---------
query = '"generative ai" AND organization AND transformation'
year_start = 2024
year_end = 2026


# Output file (you can change name and path)
output_excel_file = "C:/PYTHON/.data/openalex_sortie.xlsx"

max_results = 500
per_page = 200

config.email = "dev@atawao.com"
config.retry_count = 5
config.timeout = 30

# Ajoute cette variable dans la zone SETTINGS (comme output_excel_file)
pdf_download_dir = "C:/PYTHON/.data/openalexpdfs"   # <-- dossier où stocker les PDFs
download_pdfs = True                # mets False si tu veux désactiver

# ----------------------------------------
def safe_filename(name: str, max_len: int = 150) -> str:
    """
    Nettoie un nom de fichier pour Windows/Mac/Linux.
    """
    log("Start function: safe_filename")
    name = name.strip()
    name = re.sub(r'[<>:"/\\|?*\n\r\t]', "_", name)  # caractères interdits Windows
    name = re.sub(r"\s+", " ", name).strip()
    if len(name) > max_len:
        name = name[:max_len].rstrip()
    log("End function: safe_filename")
    return name or "file"

def guess_pdf_filename(pdf_url: str, fallback_stem: str) -> str:
    """
    Déduit un nom de fichier PDF.
    Retourne TOUJOURS une string valide.
    """
    log("Start function: guess_pdf_filename")

    try:
        parsed = urlparse(pdf_url)
        base = os.path.basename(parsed.path).strip()

        if base.lower().endswith(".pdf") and len(base) > 4:
            filename = safe_filename(base)
            log("End function: guess_pdf_filename (from url)")
            return filename

    except Exception as e:
        log(f"Filename parsing error: {e}")

    # fallback garanti
    filename = safe_filename(fallback_stem) + ".pdf"
    log("End function: guess_pdf_filename (fallback)")
    return filename


def download_pdf(pdf_url: str, out_path: Path, timeout: int = 60) -> bool:
    """
    Télécharge un PDF à partir d’un lien direct.
    Retourne True si succès, False sinon.
    """
    log("Start function: download_pdf")
    log(f"Downloading: {pdf_url}")

    try:
        headers = {
            # certains serveurs refusent sans User-Agent
            "User-Agent": "Mozilla/5.0 (PDF downloader; OpenAlex script)"
        }
        with requests.get(pdf_url, stream=True, timeout=timeout, headers=headers, allow_redirects=True) as r:
            r.raise_for_status()

            content_type = (r.headers.get("Content-Type") or "").lower()
            # ce check évite de sauvegarder une page HTML en .pdf
            if "pdf" not in content_type and not str(r.url).lower().endswith(".pdf"):
                log(f"Not a PDF (Content-Type={content_type}) -> skipped")
                log("End function: download_pdf")
                return False

            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 128):
                    if chunk:
                        f.write(chunk)

        log(f"Saved: {out_path}")
        log("End function: download_pdf")
        return True

    except Exception as e:
        log(f"Download failed: {e}")
        log("End function: download_pdf")
        return False

    filename = safe_filename(fallback_stem) + ".pdf"
    log("End function: guess_pdf_filename (fallback)")
    return filename

def download_available_pdfs(df: pd.DataFrame, pdf_dir: str) -> None:
    """
    Télécharge tous les PDFs disponibles en utilisant la colonne df['pdf_url'].
    Le dossier de sortie est donné par la variable pdf_dir.
    """
    log("Start function: download_available_pdfs")

    out_dir = Path(pdf_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if "pdf_url" not in df.columns:
        log("Column 'pdf_url' not found -> nothing to download")
        log("End function: download_available_pdfs")
        return

    total = len(df)
    count_has_link = int(df["pdf_url"].notna().sum())
    log(f"Rows: {total}, rows with pdf_url: {count_has_link}")

    ok = 0
    skipped = 0
    failed = 0

    for i, row in df.iterrows():
        pdf_url = row.get("pdf_url")
        if not pdf_url or not isinstance(pdf_url, str) or not pdf_url.strip():
            skipped += 1
            continue

        title = row.get("title") or "untitled"
        year = row.get("publication_year") or ""
        doi = row.get("doi") or ""
        fallback_stem = f"{year} - {title}"
        if doi:
            fallback_stem += f" - {doi}"

        filename = guess_pdf_filename(pdf_url, fallback_stem)
        out_path = out_dir / filename

        # évite de re-télécharger si déjà présent
        if out_path.exists() and out_path.stat().st_size > 0:
            log(f"Already exists -> skip: {out_path.name}")
            ok += 1
            continue

        success = download_pdf(pdf_url, out_path)
        if success:
            ok += 1
        else:
            failed += 1

    log(f"PDF download summary: ok={ok}, skipped_no_link={skipped}, failed={failed}")
    log("End function: download_available_pdfs")

def inverted_index_to_text(
    abstract_inverted_index: Optional[Dict[str, List[int]]]
) -> Optional[str]:
    log("Start function: inverted_index_to_text")

    if not abstract_inverted_index:
        log("No abstract_inverted_index found")
        return None

    max_pos = -1
    for positions in abstract_inverted_index.values():
        if positions:
            max_pos = max(max_pos, max(positions))

    if max_pos < 0:
        log("Abstract index empty after scan")
        return None

    words = [""] * (max_pos + 1)
    for word, positions in abstract_inverted_index.items():
        for p in positions:
            if 0 <= p < len(words):
                words[p] = word

    text = " ".join(w for w in words if w).strip()
    log("End function: inverted_index_to_text")

    return text or None


def looks_like_pdf(url: str) -> bool:
    log("Start function: looks_like_pdf")

    u = url.lower()
    result = u.endswith(".pdf") or bool(re.search(r"\.pdf(\?|$)", u))

    log(f"End function: looks_like_pdf (result={result})")
    return result


def extract_pdf_url(work: Dict[str, Any]) -> Optional[str]:
    log("Start function: extract_pdf_url")

    primary = work.get("primary_location") or {}
    pdf = primary.get("pdf_url")
    if pdf:
        log("PDF found in primary_location")
        log("End function: extract_pdf_url")
        return pdf

    for loc in (work.get("locations") or []):
        pdf = (loc or {}).get("pdf_url")
        if pdf:
            log("PDF found in locations[]")
            log("End function: extract_pdf_url")
            return pdf

    oa_url = (work.get("open_access") or {}).get("oa_url")
    if oa_url and looks_like_pdf(oa_url):
        log("PDF inferred from open_access.oa_url")
        log("End function: extract_pdf_url")
        return oa_url

    log("No PDF URL found")
    log("End function: extract_pdf_url")
    return None


def extract_landing_page_url(work: Dict[str, Any]) -> Optional[str]:
    log("Start function: extract_landing_page_url")

    primary = work.get("primary_location") or {}
    if primary.get("landing_page_url"):
        log("Landing page found in primary_location")
        log("End function: extract_landing_page_url")
        return primary["landing_page_url"]

    for loc in (work.get("locations") or []):
        if (loc or {}).get("landing_page_url"):
            log("Landing page found in locations[]")
            log("End function: extract_landing_page_url")
            return loc["landing_page_url"]

    log("No landing page URL found")
    log("End function: extract_landing_page_url")
    return None


def fetch_works(query: str, max_results: int, per_page: int, year_start: int, year_end: int):
    log("Start function: fetch_works")
    log(f"Query='{query}', max_results={max_results}, per_page={per_page}, years={year_start}-{year_end}")

    results = []

    from_date = f"{year_start}-01-01"
    to_date = f"{year_end}-12-31"

    # paginate(per_page=...) est la méthode recommandée dans pyalex
    pager = (
        Works()
        .search(query)
        .filter(
            from_publication_date=from_date,
            to_publication_date=to_date
        )
        .paginate(per_page=per_page)
    )

    for page in pager:
        log(f"Fetched a page with {len(page)} works")
        for work in page:
            results.append(work)
            if len(results) >= max_results:
                log(f"Reached max_results={max_results}")
                log("End function: fetch_works")
                return results

    log(f"Fetched {len(results)} works")
    log("End function: fetch_works")
    return results

def simplify_for_excel(work: Dict[str, Any]) -> Dict[str, Any]:
    log("Start function: simplify_for_excel")

    authors = []
    for a in (work.get("authorships") or []):
        name = ((a or {}).get("author") or {}).get("display_name")
        if name:
            authors.append(name)

    primary_location = work.get("primary_location") or {}
    source = ((primary_location.get("source") or {}) or {}).get("display_name")

    abstract_text = inverted_index_to_text(work.get("abstract_inverted_index"))

    oa = work.get("open_access") or {}

    row = {
        "id": work.get("id"),
        "doi": work.get("doi"),
        "title": work.get("title"),
        "abstract": abstract_text,
        "publication_year": work.get("publication_year"),
        "type": work.get("type"),
        "language": work.get("language"),
        "source": source,
        "pdf_url": extract_pdf_url(work),
        "landing_page_url": extract_landing_page_url(work),
        "oa_status": oa.get("oa_status"),
        "oa_url": oa.get("oa_url"),
        "cited_by_count": work.get("cited_by_count"),
        "authors": "; ".join(authors) if authors else None,
    }

    log("End function: simplify_for_excel")
    return row


def main() -> None:
    log("Lancement du programme")

    raw = fetch_works(
        query=query,
        max_results=max_results,
        per_page=per_page,
        year_start=year_start,
        year_end=year_end,
    )

    log("Transforming raw works into flat rows")

    rows = [simplify_for_excel(w) for w in raw]
    df = pd.DataFrame(rows)

    output_path = Path(output_excel_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_excel(output_path, index=False)
    if download_pdfs:
        download_available_pdfs(df, pdf_download_dir)

    log(f"Excel file saved to {output_path}")
    log(f"Total rows written: {len(df)}")

    log("Fin du programme")


if __name__ == "__main__":
    main()
