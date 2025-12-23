# -*- coding: utf-8 -*-
import os
import re
import time
import json
import shutil
import string
import requests
import pandas as pd

from bs4 import BeautifulSoup
from habanero import Crossref
from io import BytesIO
from slugify import slugify
from pdfminer.high_level import extract_text as pdf_extract_text

try:
    import trafilatura
    HAS_TRAF = True
except Exception:
    HAS_TRAF = False

# ---------------------------
# Configuration utilisateur
# ---------------------------
UNPAYWALL_EMAIL = "dev@atawao.com"  # <-- Mets ton email ici
CR = Crossref(mailto="dev@atawao.com")  # mets ton email
INPUT_FILE = r"C:/DATA/code/.params/listeDOI.xlsx"  # Excel/CSV avec une colonne 'DOI'
INPUT_COL = "DOI"                                    # nom de la colonne DOI
FALLBACK_DOIS = [
    # "10.1038/s41586-019-1666-5",
]
OUTPUT_DIR = "C:/DATA/code/.data/fulltext"
RATE_LIMIT_SECONDS = 0.5  # politesse API

UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/126.0 Safari/537.36"
)

# ---------------------------
# Helpers
# ---------------------------

def normalize_doi(doi: str) -> str | None:
    if not doi:
        return None
    s = str(doi).strip()
    if not s or s in {"?", "nan", "None"}:
        return None
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.I)
    s = re.sub(r"^doi:\s*", "", s, flags=re.I)
    return s

def safe_filename_from(meta: dict, doi: str) -> str:
    parts = []
    year = meta.get("year") or meta.get("published_date") or ""
    if isinstance(year, str):
        year = re.findall(r"\d{4}", year)
        year = year[0] if year else ""
    elif isinstance(year, int):
        year = str(year)

    title = meta.get("title") or ""
    journal = meta.get("journal_name") or meta.get("journal") or ""

    if year: parts.append(year)
    if title: parts.append(title)
    elif journal: parts.append(journal)
    parts.append(doi)

    base = " - ".join(parts)
    base = slugify(base, max_length=140, lowercase=True)
    if not base:
        base = slugify(doi)
    return base + ".txt"

def get_unpaywall_json(doi: str) -> dict | None:
    url = f"https://api.unpaywall.org/v2/{requests.utils.quote(doi)}"
    params = {"email": UNPAYWALL_EMAIL}
    r = requests.get(url, params=params, timeout=20)
    if r.status_code == 200:
        return r.json()
    return None

def pick_best_pdf_url(ujson: dict) -> tuple[str | None, dict]:
    if not ujson:
        return None, {}
    meta = {
        "title": ujson.get("title"),
        "year": ujson.get("year"),
        "journal_name": ujson.get("journal_name"),
        "is_oa": ujson.get("is_oa"),
        "oa_status": ujson.get("oa_status"),
        "license": ujson.get("license"),
    }
    if not ujson.get("is_oa"):
        return None, meta
    best = ujson.get("best_oa_location") or {}
    pdf_url = best.get("url_for_pdf")
    if pdf_url:
        return pdf_url, meta
    for loc in (ujson.get("oa_locations") or []):
        if loc.get("url_for_pdf"):
            return loc["url_for_pdf"], meta
    return None, meta

def pick_best_html_url(ujson: dict) -> str | None:
    if not ujson or not ujson.get("is_oa"):
        return None
    best = ujson.get("best_oa_location") or {}
    if best.get("url"):
        return best["url"]
    for loc in (ujson.get("oa_locations") or []):
        if loc.get("url"):
            return loc["url"]
    return None

def download_pdf_text(pdf_url: str) -> str | None:
    """Télécharge un PDF et extrait le texte avec pdfminer.six."""
    try:
        headers = {"User-Agent": UA, "Referer": "https://doi.org/"}
        r = requests.get(pdf_url, headers=headers, timeout=60, stream=True, allow_redirects=True)
        if r.status_code >= 400:
            return None
        content = r.content  # parfois le header n'indique pas application/pdf
        text = pdf_extract_text(BytesIO(content))
        text = re.sub(r"\n{3,}", "\n\n", text or "")
        return text.strip() if text else None
    except Exception:
        return None

def download_html_text(page_url: str) -> str | None:
    if not HAS_TRAF:
        return None
    try:
        downloaded = trafilatura.fetch_url(page_url, timeout=60)
        if not downloaded:
            return None
        text = trafilatura.extract(downloaded, include_comments=False, include_tables=False)
        return text.strip() if text else None
    except Exception:
        return None

def download_html_text_with_headers(page_url: str) -> str | None:
    """Version HTML avec headers pour contourner consent pages/cdn."""
    if not HAS_TRAF:
        return None
    try:
        headers = {"User-Agent": UA, "Referer": "https://doi.org/", "Accept-Language": "en,fr;q=0.9"}
        r = requests.get(page_url, headers=headers, timeout=60)
        if r.status_code >= 400:
            return None
        html = r.text
        text = trafilatura.extract(html, include_comments=False, include_tables=False, url=page_url)
        return text.strip() if text else None
    except Exception:
        return None

def extract_source_html_text_as_last_resort(doi: str, timeout=60) -> tuple[str | None, str | None]:
    """
    Dernier recours : récupère le texte depuis la page HTML de la source éditeur.
    1) Résout https://doi.org/<doi> → landing
    2) Tente trafilatura avec headers
    3) Si échec, fallback BeautifulSoup sur blocs d'article connus
    Retourne (texte, landing_url) ou (None, None)
    """
    landing = resolve_doi_landing(doi)
    if not landing:
        return None, None

    # 1) Trafilatura avec headers (plus robuste que fetch_url)
    text = download_html_text_with_headers(landing)
    if text:
        return text, landing

    # 2) Fallback BeautifulSoup : extraire les blocs de contenu
    try:
        headers = {"User-Agent": UA, "Referer": "https://doi.org/", "Accept-Language": "en,fr;q=0.9"}
        r = requests.get(landing, headers=headers, timeout=timeout)
        if r.status_code >= 400:
            return None, landing
        soup = BeautifulSoup(r.text, "lxml")

        # Sélecteurs fréquents pour corps d’article (Nature/Springer/Wiley/Elsevier/T&F/Frontiers)
        candidates = [
            # Nature / Springer
            "div.c-article-body", "article.c-article-body", "div#Abs1-content",
            # Wiley
            "div.article-section__content", "div.article-section", "section#main-content",
            # Elsevier / ScienceDirect
            "div#body", "div.article__body", "div#main_content", "div#frag_1", "div.abstract", "div#abs",
            # Taylor & Francis
            "div#content", "div.NLM_sec", "div.article",
            # fallback génériques
            "article", "main", "div#content", "div.content",
        ]

        # Concatène le texte des meilleurs candidats
        collected = []
        seen = set()
        for sel in candidates:
            for node in soup.select(sel):
                txt = node.get_text(separator="\n", strip=True)
                if txt and txt not in seen:
                    seen.add(txt)
                    collected.append(txt)

        # Si rien via sélecteurs, fallback global (risque de bruit, on filtre un peu)
        if not collected:
            body = soup.find("body")
            if body:
                txt = body.get_text(separator="\n", strip=True)
                # petit nettoyage : réduire les sauts de ligne multiples
                txt = re.sub(r"\n{3,}", "\n\n", txt or "")
                if txt:
                    collected.append(txt)

        if collected:
            # garde le plus long bloc (souvent le corps)
            collected.sort(key=len, reverse=True)
            return collected[0], landing

    except Exception:
        pass

    return None, landing

def write_txt(path: str, header: dict, body: str):
    os.makedirs(os.path.dirname(path), exist_ok=True
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write("### METADATA ###\n")
        for k, v in header.items():
            f.write(f"{k}: {v}\n")
        f.write("\n### FULLTEXT ###\n\n")
        f.write(body)

def resolve_doi_landing(doi: str, timeout=20) -> str | None:
    """Equivalent du formulaire doi.org : retourne l'URL de la page éditeur."""
    try:
        r = requests.get(
            f"https://doi.org/{doi}",
            headers={"User-Agent": UA},
            allow_redirects=True,
            timeout=timeout,
        )
        if 200 <= r.status_code < 400:
            return r.url
    except Exception:
        pass
    return None

def resolve_doi_pdf_via_content_negotiation(doi: str, timeout=20) -> str | None:
    """Tente d'obtenir directement un PDF via content negotiation."""
    try:
        r = requests.get(
            f"https://doi.org/{doi}",
            headers={"User-Agent": UA, "Accept": "application/pdf"},
            allow_redirects=True,
            timeout=timeout,
        )
        if r.url.lower().endswith(".pdf") or "application/pdf" in (r.headers.get("Content-Type") or "").lower():
            return r.url
    except Exception:
        pass
    return None

def resolve_doi_crossref_json(doi: str, timeout=20) -> dict | None:
    """Récupère les métadonnées Crossref via content negotiation."""
    try:
        r = requests.get(
            f"https://doi.org/{doi}",
            headers={"User-Agent": UA, "Accept": "application/vnd.crossref+json"},
            timeout=timeout,
        )
        if r.status_code == 200:
            return r.json()
    except Exception:
        pass
    return None

def resolve_pdf_via_crossref(doi: str) -> str | None:
    """Cherche un lien PDF exposé par l'éditeur via Crossref (habanero)."""
    try:
        r = CR.works(ids=doi)
        msg = r.get("message", {})
        for link in msg.get("link", []):
            ct = (link.get("content-type") or "").lower()
            url = link.get("URL")
            if url and "pdf" in ct:
                return url
    except Exception:
        pass
    return None

def find_pdf_on_publisher_page(page_url: str, timeout=30) -> str | None:
    """Trouve un lien PDF en parsant la page éditeur (meta/link/a .pdf)."""
    try:
        headers = {"User-Agent": UA, "Referer": "https://doi.org/", "Accept-Language": "en,fr;q=0.9"}
        r = requests.get(page_url, headers=headers, timeout=timeout)
        if r.status_code >= 400:
            return None
        soup = BeautifulSoup(r.text, "lxml")
        meta_pdf = soup.find("meta", attrs={"name": "citation_pdf_url"})
        if meta_pdf and meta_pdf.get("content"):
            return meta_pdf["content"].strip()
        for a in soup.find_all("a", href=True):
            href = a["href"]
            if ".pdf" in href.lower():
                return requests.compat.urljoin(page_url, href)
        for l in soup.find_all("link", href=True):
            href = l["href"]
            if ".pdf" in href.lower():
                return requests.compat.urljoin(page_url, href)
    except Exception:
        return None
    return None

def sciencedirect_pdf_from_landing(landing_url: str, timeout=40) -> str | None:
    """
    Cas particulier ScienceDirect (Elsevier) : construit et teste des URLs PDF
    à partir du PII et/ou des liens présents dans la page.
    """
    if "sciencedirect.com" not in (landing_url or ""):
        return None

    headers = {"User-Agent": UA, "Referer": landing_url, "Accept-Language": "en,fr;q=0.9"}

    # 1) Extraire le PII de l'URL
    m = re.search(r"/pii/([A-Z0-9]+)", landing_url)
    pii = m.group(1) if m else None
    candidates = []

    # 2) Générer les patterns PDF connus
    if pii:
        base = f"https://www.sciencedirect.com/science/article/pii/{pii}"
        candidates += [
            base + "/pdf",
            base + "/pdf?download=true",
            base + "/pdfft",
            base + "/pdfft?download=true",
            base + "/pdfft?isDTMRedir=true&download=true",
        ]

    # 3) Scraper la page pour récolter d'autres liens .pdf
    try:
        r = requests.get(landing_url, headers=headers, timeout=timeout)
        if r.status_code < 400:
            soup = BeautifulSoup(r.text, "lxml")
            for a in soup.find_all("a", href=True):
                href = a["href"]
                if ".pdf" in href.lower():
                    candidates.append(requests.compat.urljoin(landing_url, href))
            for l in soup.find_all("link", href=True):
                href = l["href"]
                if ".pdf" in href.lower():
                    candidates.append(requests.compat.urljoin(landing_url, href))
    except Exception:
        pass

    # 4) Tester les candidates
    seen = set()
    for url in candidates:
        if not url or url in seen:
            continue
        seen.add(url)
        try:
            h = requests.head(url, headers=headers, allow_redirects=True, timeout=timeout)
            ctype = (h.headers.get("Content-Type") or "").lower()
            if "pdf" not in ctype:
                g = requests.get(url, headers=headers, stream=True, allow_redirects=True, timeout=timeout)
                ctype = (g.headers.get("Content-Type") or "").lower()
                g.close()
            if "pdf" in ctype or url.lower().endswith(".pdf"):
                return url
        except Exception:
            continue

    return None

# ---------------------------
# Entrée : collecte des DOIs
# ---------------------------
def load_dois() -> list[str]:
    dois = set()
    if os.path.exists(INPUT_FILE):
        ext = os.path.splitext(INPUT_FILE)[1].lower()
        if ext in {".xlsx", ".xls"}:
            df = pd.read_excel(INPUT_FILE)  # nécessite openpyxl pour .xlsx
        else:
            df = pd.read_csv(INPUT_FILE)
        cols = {c.strip().lower(): c for c in df.columns}
        if INPUT_COL.lower() not in cols:
            if df.shape[1] == 1:
                colname = df.columns[0]
            else:
                raise ValueError(f"Colonne '{INPUT_COL}' introuvable dans {INPUT_FILE}. Colonnes: {list(df.columns)}")
        else:
            colname = cols[INPUT_COL.lower()]
        for v in df[colname].dropna().tolist():
            doi = normalize_doi(v)
            if doi:
                dois.add(doi)
    for v in FALLBACK_DOIS:
        doi = normalize_doi(v)
        if doi:
            dois.add(doi)
    return sorted(dois)

# ---------------------------
# Boucle principale
# ---------------------------
def main():
    if not UNPAYWALL_EMAIL or "@" not in UNPAYWALL_EMAIL:
        raise SystemExit("⚠️ Renseigne UNPAYWALL_EMAIL (ex: ton.email@domaine.com).")

    dois = load_dois()
    if not dois:
        print("Aucun DOI fourni (INPUT_FILE vide et FALLBACK_DOIS vide).")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"DOIs à traiter: {len(dois)}")

    for idx, doi in enumerate(dois, start=1):
        print(f"[{idx}/{len(dois)}] {doi} …")
        try:
            # 1) Unpaywall
            uj = get_unpaywall_json(doi)
            if not uj or not uj.get("is_oa"):
                print("   - Pas en libre accès (Unpaywall).")
                time.sleep(RATE_LIMIT_SECONDS)
                continue

            pdf_url, meta = pick_best_pdf_url(uj)
            html_url = pick_best_html_url(uj)

            # 2) Récupération primaire (Unpaywall)
            body_text = None
            src_url = None
            if pdf_url:
                body_text = download_pdf_text(pdf_url)
                src_url = pdf_url
                if body_text:
                    print("   - PDF OA récupéré.")
            if not body_text and html_url:
                body_text = download_html_text(html_url)
                src_url = html_url
                if body_text:
                    print("   - HTML OA récupéré.")

            # 3) Fallbacks éditeur (content negotiation → Crossref → landing → cas ScienceDirect → HTML avec headers)
            if not body_text:
                # 3.0 doi.org -> content negotiation pour PDF direct
                pdf_url_cn = resolve_doi_pdf_via_content_negotiation(doi)
                if pdf_url_cn:
                    txt = download_pdf_text(pdf_url_cn)
                    if txt:
                        body_text = txt
                        src_url = pdf_url_cn
                        print("   - PDF via doi.org (content negotiation) récupéré.")

            if not body_text:
                # 3.1 Crossref : lien PDF du publisher (si dispo)
                pdf_url2 = resolve_pdf_via_crossref(doi)
                if pdf_url2:
                    txt = download_pdf_text(pdf_url2)
                    if txt:
                        body_text = txt
                        src_url = pdf_url2
                        print("   - PDF éditeur (Crossref) récupéré.")

            if not body_text:
                # 3.2 doi.org → page éditeur → chercher PDF
                landing = resolve_doi_landing(doi)
                if landing:
                    # 3.2.a ScienceDirect (Elsevier) : essais spécifiques
                    if "sciencedirect.com" in landing:
                        sd_pdf = sciencedirect_pdf_from_landing(landing)
                        if sd_pdf:
                            txt = download_pdf_text(sd_pdf)
                            if txt:
                                body_text = txt
                                src_url = sd_pdf
                                print("   - PDF éditeur (ScienceDirect) récupéré.")

                    # 3.2.b Générique : meta/link .pdf
                    if not body_text:
                        pdf_url3 = find_pdf_on_publisher_page(landing)
                        if pdf_url3:
                            txt = download_pdf_text(pdf_url3)
                            if txt:
                                body_text = txt
                                src_url = pdf_url3
                                print("   - PDF éditeur (via page d'atterrissage) récupéré.")

                    # 3.3 Fallback HTML avec headers
                    if not body_text and HAS_TRAF:
                        html_text = download_html_text_with_headers(landing)
                        if html_text:
                            body_text = html_text
                            src_url = landing
                            print("   - Texte HTML éditeur récupéré (headers).")

                    # 3.4 Dernier recours : extraire le HTML de la source coûte que coûte
                    if not body_text:
                        last_txt, last_src = extract_source_html_text_as_last_resort(doi)
                        if last_txt:
                            body_text = last_txt
                            src_url = last_src
                            print("   - Texte HTML éditeur récupéré (last resort).")

            time.sleep(0.5)  # politesse

            if not body_text:
                print("   - Échec extraction texte (après fallbacks).")
                time.sleep(RATE_LIMIT_SECONDS)
                continue

            # 4) Écriture fichier
            header = {
                "doi": doi,
                "source_url": src_url or "",
                "title": meta.get("title") or "",
                "year": meta.get("year") or "",
                "journal": meta.get("journal_name") or "",
                "oa_status": meta.get("oa_status") or "",
                "license": meta.get("license") or "",
            }
            fname = safe_filename_from(meta, doi)
            outpath = os.path.join(OUTPUT_DIR, fname)
            write_txt(outpath, header, body_text)
            print(f"   → OK: {outpath}")

            time.sleep(RATE_LIMIT_SECONDS)

        except KeyboardInterrupt:
            print("\nInterrompu par l'utilisateur.")
            break
        except Exception as e:
            print(f"   - Erreur: {e}")
            time.sleep(RATE_LIMIT_SECONDS)
            continue

if __name__ == "__main__":
    main()

