# -*- coding: utf-8 -*-

# ========================
# Bibliothèques
# ========================
import os
import time
import pandas as pd
from datetime import datetime, timedelta

from Bio import Entrez, Medline
from datetime import datetime
from colorama import Fore

import re
from xml.etree import ElementTree as ET

try:
    import metapub
    HAS_METAPUB = True
except Exception:
    HAS_METAPUB = False


# ========================
# Paramètres
# ========================
BATCH_SIZE = 50  # taille de lot pour efetch (<=100 recommandé)
TERME_CHERCHE = "(('Colorectal Neoplasms'[MeSH Terms] OR colorectal cancer[Title/Abstract] OR colon cancer[Title/Abstract] OR rectal cancer[Title/Abstract]) AND ('Artificial Intelligence'[MeSH Terms] OR artificial intelligence[Title/Abstract] OR machine learning[Title/Abstract] OR deep learning[Title/Abstract] OR neural network*[Title/Abstract]))"
DATE_DEBUT = "2025/09/01"
DATE_FIN = "2026/03/01"
MAIL_PUBMED = "dev@atawao.com"
BASE_NCBI = "pubmed"
ORDRE_TRI = "pub date"
FORMAT_RESULTAT = "medline"
FORMAT_ENREGISTREMENT = "text"
MAX_ESearch_RESULTS = 9000  # max de résultats par sous-requête (sécurité < 9999)



APIKEY = "2bad5e3552250d73f6000974773fb3d05c08"
os.environ["NCBI_API_KEY"] = APIKEY

fichier_facteur_impact = "C:/PYTHON/.params/impact_factor.txt"
fichier_pubmed_sortie = "C:/PYTHON/.data/pubmed_sortie.xlsx"

# ========================
# Fonctions utilitaires
# ========================
def ensure_pmc_prefix(pmcid: str) -> str:
    if not pmcid:
        return ""
    pmcid = str(pmcid).strip()
    return pmcid if pmcid.upper().startswith("PMC") else f"PMC{pmcid}"

def fetch_pmc_fulltext_xml(pmcid: str) -> str:
    """Retourne le XML JATS complet depuis PMC (si disponible)."""
    pmcid = ensure_pmc_prefix(pmcid)
    if not pmcid:
        return ""

    try:
        h = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
        xml = h.read()
        time.sleep(0.34)  # rester gentil avec NCBI
        return xml
    except Exception as e:
        print(Fore.YELLOW + f"[WARN] Impossible de récupérer le full text PMC {pmcid}: {e}" + Fore.RESET)
        return ""

def jats_xml_to_plain_text(xml_str: str, max_chars: int = 200000) -> str:
    """
    Extrait un texte lisible depuis du XML JATS.
    - On essaye de prendre surtout <body>.
    - On limite la taille pour éviter un Excel énorme.
    """
    if not xml_str:
        return ""

    try:
        root = ET.fromstring(xml_str)
    except Exception:
        # XML invalide / réponse inattendue
        return ""

    # Cherche un <body> (souvent dans JATS)
    body = None
    for el in root.iter():
        if el.tag.lower().endswith("body"):
            body = el
            break

    target = body if body is not None else root

    # Concatène tous les textes
    text = " ".join(t.strip() for t in target.itertext() if t and t.strip())
    text = re.sub(r"\s+", " ", text).strip()

    if max_chars and len(text) > max_chars:
        text = text[:max_chars] + " …[TRONQUÉ]"
    return text

def now():
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")

def parse_yyyymmdd(s: str) -> datetime:
    return datetime.strptime(s, "%Y/%m/%d")

def fmt_yyyymmdd(d: datetime) -> str:
    return d.strftime("%Y/%m/%d")

def build_term(base_term: str, d1: datetime, d2: datetime) -> str:
    # PDAT inclusif sur bornes en pratique (NCBI accepte bien ce format)
    return f'{base_term} AND ({fmt_yyyymmdd(d1)}[PDAT] : {fmt_yyyymmdd(d2)}[PDAT])'

def esearch_count(term: str) -> int:
    h = Entrez.esearch(db=BASE_NCBI, term=term, sort=ORDRE_TRI, retmax=1)
    r = Entrez.read(h)
    return int(r.get("Count", 0))

def fetch_all_ids_for_term(term: str) -> list[str]:
    ids = []
    # Ici retstart reste toujours < 9999 car on garantit Count <= MAX_ESearch_RESULTS
    retstart = 0
    while True:
        h = Entrez.esearch(
            db=BASE_NCBI,
            term=term,
            sort=ORDRE_TRI,
            retstart=retstart,
            retmax=BATCH_SIZE
        )
        r = Entrez.read(h)
        batch = r.get("IdList", [])
        if not batch:
            break
        ids.extend(batch)
        retstart += len(batch)
        time.sleep(0.34)
    return ids

def split_date_ranges(base_term: str, d1: datetime, d2: datetime) -> list[tuple[datetime, datetime]]:
    """
    Découpe [d1, d2] en sous-plages de dates pour que chaque sous-requête ait <= MAX_ESearch_RESULTS.
    Méthode : on coupe en 2 (bisection) jusqu'à passer sous le seuil.
    """
    term = build_term(base_term, d1, d2)
    cnt = esearch_count(term)

    if cnt == 0:
        return []

    if cnt <= MAX_ESearch_RESULTS:
        return [(d1, d2)]

    # Si on est déjà au jour près et c'est encore trop, on ne peut pas couper plus fin
    if d1.date() == d2.date():
        # dans ce cas, il faudra une autre stratégie (ex: ajouter filtre supplémentaire), mais c’est rare
        return [(d1, d2)]

    mid = d1 + (d2 - d1) / 2
    mid = datetime(mid.year, mid.month, mid.day)  # on tronque à minuit pour éviter les demi-jours

    left = split_date_ranges(base_term, d1, mid)
    right = split_date_ranges(base_term, mid + timedelta(days=1), d2)
    return left + right

def normalize_journal_name(name: str) -> str:
    if not isinstance(name, str):
        return ""
    return name.strip().lower()


def extract_ids_from_xml_records(xml_records):
    mapping = {}
    for art in xml_records:
        try:
            pmid = str(art["MedlineCitation"]["PMID"])
            doi = None
            pmcid = None
            other_ids = []
            for aid in art.get("PubmedData", {}).get("ArticleIdList", []):
                id_val = str(aid)
                id_type = getattr(aid, "attributes", {}).get("IdType", "").lower()
                if id_type == "doi":
                    doi = id_val
                elif id_type in ("pmcid", "pmc"):
                    pmcid = id_val
                else:
                    other_ids.append(id_val)
            mapping[pmid] = {"doi": doi, "pmcid": pmcid, "other_ids": other_ids}
        except Exception:
            continue
    return mapping


def safe_metapub_doi(pmid: str) -> str | None:
    if not HAS_METAPUB:
        return None
    try:
        fi = metapub.FindIt(pmid=pmid)
        return getattr(fi, "doi", None)
    except Exception as e:
        print(Fore.YELLOW + f"[WARN] FindIt failed for PMID {pmid}: {e}" + Fore.RESET)
        return None


def list_to_str(x):
    if x == "?":
        return "?"
    if isinstance(x, list):
        return ", ".join(map(str, x))
    return str(x) if x is not None else ""


# ========================
# Initialisation Entrez
# ========================
Entrez.email = MAIL_PUBMED
Entrez.api_key = os.environ["NCBI_API_KEY"]

# ========================
# Lecture facteur d'impact (robuste)
# ========================
print(f"{now()} - Ouverture du fichier des facteurs d'impact")

def _clean_col(s: str) -> str:
    # normalise les noms de colonnes (enlève BOM, espaces, accents basiques)
    import re, unicodedata
    s = (s or "").replace("\ufeff", "")  # enlève BOM éventuel
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    s = s.strip().lower()
    s = re.sub(r"\s+", "_", s)
    return s

# Essais de lecture robustes
df_impact = pd.read_csv(
    fichier_facteur_impact,
    sep=";",
    header=0,
    encoding="utf-8-sig",   # retire un éventuel BOM
    engine="python",
)

# Nettoie les noms de colonnes
df_impact.columns = [_clean_col(c) for c in df_impact.columns]

# Si le fichier ne contient pas d'entête correcte, on tente en header=None
if not set(df_impact.columns) - {0, 1} and len(df_impact.columns) <= 3:
    df_impact = pd.read_csv(
        fichier_facteur_impact,
        sep=";",
        header=None,
        encoding="utf-8-sig",
        engine="python",
    )
    df_impact.columns = [_clean_col(str(c)) for c in df_impact.columns]

# Détecte la colonne "titre" (journal) et la colonne "if" (impact factor) via synonymes
titre_candidates = {"titre", "journal", "revue", "title", "journal_title", "source_title", "name"}
if_candidates   = {"if", "impact_factor", "impactfactor", "impact", "jif"}

def _pick(colset, candidates):
    for cand in candidates:
        if cand in colset:
            return cand
    # tentative de recherche floue simple
    for c in colset:
        for cand in candidates:
            if cand in c:
                return c
    return None

colset = set(df_impact.columns)

col_titre = _pick(colset, titre_candidates)
col_if    = _pick(colset, if_candidates)

# Cas particulier: si seulement 2 colonnes, on suppose [titre, if]
if col_titre is None or col_if is None:
    if len(df_impact.columns) == 2:
        df_impact.columns = ["titre", "if"]
        col_titre, col_if = "titre", "if"
    elif len(df_impact.columns) == 3:
        # certains fichiers ont un index en 1re col.
        df_impact.columns = ["index", "titre", "if"]
        col_titre, col_if = "titre", "if"

# Si toujours introuvable -> message explicite
if col_titre is None or col_if is None:
    raise ValueError(
        f"Impossible d'identifier les colonnes 'titre' (journal) et 'if' dans {fichier_facteur_impact}. "
        f"Colonnes trouvées: {list(df_impact.columns)}\n"
        f"Attendu (synonymes acceptés): titre/journal/revue/title et if/impact_factor."
    )

# Normalisation pour le join avec Medline
df_impact["titre"] = df_impact[col_titre].astype(str)
df_impact["if"]    = pd.to_numeric(df_impact[col_if], errors="coerce")
df_impact["titre_norm"] = df_impact["titre"].str.strip().str.lower()

# ========================
# Recherche des PMIDs (tous)
# ========================
BASE_TERM = TERME_CHERCHE
d1 = parse_yyyymmdd(DATE_DEBUT)
d2 = parse_yyyymmdd(DATE_FIN)

TERM = build_term(BASE_TERM, d1, d2)
print(f"{now()} - Requête PubMed: {TERM}")

# Total global (info)
total_count = esearch_count(TERM)
print(f"{now()} - Nombre total de publications : {total_count}")

if total_count == 0:
    print(Fore.YELLOW + f"{now()} - Aucun résultat trouvé. Fin." + Fore.RESET)
    raise SystemExit(0)

# 1) Découpe en sous-requêtes par dates (chaque sous-requête <= 9000)
ranges = split_date_ranges(BASE_TERM, d1, d2)
print(f"{now()} - Nombre de sous-plages de dates : {len(ranges)}")

# 2) Récupère tous les PMIDs sous-plage par sous-plage
all_ids = []
for (ra, rb) in ranges:
    sub_term = build_term(BASE_TERM, ra, rb)
    sub_count = esearch_count(sub_term)
    print(f"{now()} - Sous-requête {fmt_yyyymmdd(ra)} -> {fmt_yyyymmdd(rb)} : {sub_count} résultats")

    ids_sub = fetch_all_ids_for_term(sub_term)
    all_ids.extend(ids_sub)

    # petite pause entre sous-requêtes
    time.sleep(1)

print(f"{now()} - Total d'identifiants récupérés : {len(all_ids)}")

# ========================
# Collecte des données en lots de 50
# ========================
rows = []
for i in range(0, len(all_ids), BATCH_SIZE):
    batch_ids = all_ids[i:i+BATCH_SIZE]
    print(f"{now()} - Lecture détails {i+1}-{i+len(batch_ids)} / {len(all_ids)}")

    # MEDLINE pour texte + XML pour DOI
    h_med = Entrez.efetch(db=BASE_NCBI, id=",".join(batch_ids),
                          rettype=FORMAT_RESULTAT, retmode=FORMAT_ENREGISTREMENT)
    records = list(Medline.parse(h_med))

    time.sleep(0.5)
    h_xml = Entrez.efetch(db=BASE_NCBI, id=",".join(batch_ids), retmode="xml")
    xml_records = Entrez.read(h_xml)
    articles = xml_records.get("PubmedArticle", [])
    pmid_to_ids = extract_ids_from_xml_records(articles)


    for record in records:
        try:
            pmid = record.get("PMID", "?")
            journal_medline = record.get("JT", "?")
            journal_norm = normalize_journal_name(journal_medline)

            facteur_impact = ""
            if journal_norm:
                match = df_impact.loc[df_impact["titre_norm"] == journal_norm, "if"]
                if not match.empty:
                    facteur_impact = match.iloc[0]

            ids_bundle = pmid_to_ids.get(str(pmid), {}) if pmid else {}
            doi = ids_bundle.get("doi")
            pmcid = ids_bundle.get("pmcid")

            if not doi:
                # Cherche DOI dans LID / AID
                lid = record.get("LID", "")
                aid = record.get("AID", "")
                all_candidates = []
                if isinstance(lid, list):
                    all_candidates.extend(lid)
                elif isinstance(lid, str):
                    all_candidates.append(lid)
                if isinstance(aid, list):
                    all_candidates.extend(aid)
                elif isinstance(aid, str):
                    all_candidates.append(aid)

                for v in all_candidates:
                    if "[doi]" in str(v).lower():
                        doi = str(v).split()[0]
                        break

            if not doi and HAS_METAPUB:
                doi = safe_metapub_doi(pmid)

            auteurs = list_to_str(record.get("AU", "?"))
            titre = record.get("TI", "?")
            copyright_note = record.get("CI", "")
            date_creation = record.get("DEP") or record.get("DP") or ""
            oid = list_to_str(record.get("OID", ""))
            oci = list_to_str(record.get("OCI", ""))
            resume = record.get("AB", "")

            texte_complet = ""
            if pmcid:
                xml_full = fetch_pmc_fulltext_xml(pmcid)
                texte_complet = jats_xml_to_plain_text(xml_full)


            rows.append({
                "id": pmid,
                "titre": titre,
                "Ifactor": facteur_impact,
                "DOI": doi or "",
                "PMC": pmcid or "",
                "OID": oid,
                "OCI": oci,
                "copyright": copyright_note,
                "auteurs": auteurs,
                "journal": journal_medline or "",
                "creation": date_creation,
                "resume": resume,
                "texte_complet": texte_complet,
            })
        except Exception as e:
            print(Fore.RED + f"{now()} - Erreur PMID {record.get('PMID','?')}: {e}" + Fore.RESET)
            continue

# ========================
# Export
# ========================
print(f"{now()} - Création du DataFrame et export Excel")
df_pages = pd.DataFrame(rows, columns=[
    "id", "titre", "Ifactor", "DOI", "PMC", "OID", "OCI",
    "copyright", "auteurs", "journal", "creation", "resume", "texte_complet"
])

if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

df_pages.to_excel(fichier_pubmed_sortie, index=False)
print(f"{now()} - Terminé → {fichier_pubmed_sortie}")
