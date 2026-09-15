import requests
import pandas as pd
import time
from bs4 import BeautifulSoup
from urllib.parse import urlencode, urljoin

# ============================================================
# PARAMÈTRES
# ============================================================

NB_MAX_ABSTRACTS = 1000

# Nombre de résultats demandés par page ERS
RESULTATS_PAR_PAGE = 100

FICHIER_SORTIE = "C:/PYTHON/.data/ERS_IA_respiratoire.xlsx"

# Années à analyser
ANNEE_DEBUT = 2025
ANNEE_FIN = 2026

PAUSE_REQUETES = 0.5


# ============================================================
# TERMES IA
# ============================================================

TERMES_IA = [
    "artificial intelligence",
    "machine learning",
    "deep learning",
    "neural network",
    "neural networks",
    "large language model",
    "large language models",
    "foundation model",
    "computer vision",
    "natural language processing",
    "generative AI",
    "predictive model",
    "predictive modelling",
    "algorithm"
]


# ============================================================
# TERMES RESPIRATOIRES
# ============================================================

TERMES_RESPIRATOIRE = [
    "respiratory",
    "pulmonary",
    "lung",
    "lungs",
    "asthma",
    "copd",
    "chronic obstructive pulmonary disease",
    "bronchiectasis",
    "pulmonary fibrosis",
    "interstitial lung disease",
    "ild",
    "spirometry",
    "spirometric",
    "emphysema",
    "airway",
    "airways",
    "bronchial",
    "pneumonia",
    "sleep apnea",
    "sleep apnoea",
    "cystic fibrosis",
    "pulmonary hypertension"
]


# ============================================================
# SESSION HTTP
# ============================================================

session = requests.Session()

session.headers.update({
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 "
        "(KHTML, like Gecko) "
        "Chrome/142.0 Safari/537.36"
    ),
    "Accept-Language": "en-US,en;q=0.9"
})


BASE_URL = "https://publications.ersnet.org"

SEARCH_URL = (
    "https://publications.ersnet.org/search"
)


# ============================================================
# FONCTIONS
# ============================================================

def contient_un_terme(texte, termes):

    texte = str(texte).lower()

    return any(
        terme.lower() in texte
        for terme in termes
    )


# ============================================================
# CONSTRUCTION DE LA REQUÊTE IA
# ============================================================

REQUETE_IA = " OR ".join(
    f'"{terme}"'
    for terme in TERMES_IA
)

print("Requête IA :")
print(REQUETE_IA)


# ============================================================
# EXTRACTION DES RÉSULTATS D'UNE PAGE ERS
# ============================================================

def extraire_page(page):

    params = [
        ("query", REQUETE_IA),

        # Filtre Conference abstract
        (
            "f[0]",
            "content_type_facet:periodical-meeting-report"
        ),

        ("page", page),

        (
            "items_per_page",
            RESULTATS_PAR_PAGE
        ),

        (
            "sort_by",
            "date_ppub_facet"
        ),

        (
            "sort_order",
            "DESC"
        )
    ]


    try:

        response = session.get(
            SEARCH_URL,
            params=params,
            timeout=60
        )

        response.raise_for_status()

    except requests.RequestException as e:

        print(
            f"Erreur page {page} : {e}"
        )

        return [], False


    soup = BeautifulSoup(
        response.text,
        "html.parser"
    )


    # ========================================================
    # IDENTIFICATION DES BLOCS DE RÉSULTATS
    # ========================================================

    resultats = []

    # ERS peut faire évoluer légèrement son HTML.
    # On essaie plusieurs structures possibles.

    candidats = soup.select(
        "article, "
        ".search-result, "
        ".views-row, "
        ".node--type-periodical-meeting-report"
    )


    # Si les sélecteurs précédents ne trouvent rien,
    # on récupère les blocs contenant un titre h3/h2.
    if not candidats:

        candidats = []

        for titre in soup.select(
            "h2, h3"
        ):

            parent = titre.find_parent(
                ["article", "div", "li"]
            )

            if parent:
                candidats.append(parent)


    urls_vues = set()


    for bloc in candidats:

        texte_bloc = bloc.get_text(
            " ",
            strip=True
        )

        if not texte_bloc:
            continue


        # ----------------------------------------------------
        # Titre + URL
        # ----------------------------------------------------

        lien = bloc.select_one(
            "h2 a, h3 a, a[href*='/content/']"
        )

        if lien is None:
            continue

        titre = lien.get_text(
            " ",
            strip=True
        )

        href = lien.get(
            "href",
            ""
        )

        if not href:
            continue

        url = urljoin(
            BASE_URL,
            href
        )


        # Évite les doublons
        if url in urls_vues:
            continue

        urls_vues.add(url)


        # ----------------------------------------------------
        # Vérification Conference abstract
        # ----------------------------------------------------

        if (
            "conference abstract"
            not in texte_bloc.lower()
            and
            "suppl"
            not in texte_bloc.lower()
        ):
            continue


        resultats.append({
            "Titre": titre,
            "URL": url,
            "Texte_resultat": texte_bloc
        })


    # --------------------------------------------------------
    # Détection page suivante
    # --------------------------------------------------------

    page_suivante = False

    next_link = soup.select_one(
        'a[rel="next"], '
        '.pager__item--next a, '
        'a[title*="next" i]'
    )

    if next_link:
        page_suivante = True


    return resultats, page_suivante


# ============================================================
# EXTRACTION DE LA PAGE D'UN ABSTRACT
# ============================================================

def extraire_detail_abstract(url):

    resultat = {
        "Auteurs": "",
        "Publication": "",
        "Date_publication": "",
        "DOI": "",
        "Abstract": ""
    }


    try:

        response = session.get(
            url,
            timeout=60
        )

        response.raise_for_status()

    except requests.RequestException as e:

        print(
            f"Erreur abstract : {e}"
        )

        return resultat


    soup = BeautifulSoup(
        response.text,
        "html.parser"
    )


    # ========================================================
    # AUTEURS
    # ========================================================

    auteurs = soup.select(
        ".authors a, "
        ".author-name, "
        ".highwire-citation-authors "
        ".highwire-citation-author"
    )

    liste_auteurs = []

    for auteur in auteurs:

        texte = auteur.get_text(
            " ",
            strip=True
        )

        if texte:
            liste_auteurs.append(
                texte
            )

    resultat["Auteurs"] = " | ".join(
        dict.fromkeys(
            liste_auteurs
        )
    )


    # ========================================================
    # DOI
    # ========================================================

    doi = ""

    meta_doi = soup.find(
        "meta",
        attrs={
            "name": "citation_doi"
        }
    )

    if meta_doi:

        doi = meta_doi.get(
            "content",
            ""
        )

    resultat["DOI"] = doi


    # ========================================================
    # DATE
    # ========================================================

    meta_date = soup.find(
        "meta",
        attrs={
            "name": "citation_publication_date"
        }
    )

    if meta_date:

        resultat["Date_publication"] = (
            meta_date.get(
                "content",
                ""
            )
        )


    # ========================================================
    # PUBLICATION
    # ========================================================

    meta_journal = soup.find(
        "meta",
        attrs={
            "name": "citation_journal_title"
        }
    )

    if meta_journal:

        resultat["Publication"] = (
            meta_journal.get(
                "content",
                ""
            )
        )


    # ========================================================
    # ABSTRACT
    # ========================================================

    abstract = soup.select_one(
        ".abstract, "
        ".section.abstract, "
        "#abstract, "
        ".article-section__abstract, "
        ".highwire-markup .abstract"
    )


    if abstract:

        resultat["Abstract"] = (
            abstract.get_text(
                " ",
                strip=True
            )
        )

    else:

        # Solution de secours via meta description
        meta_abstract = soup.find(
            "meta",
            attrs={
                "name": "citation_abstract"
            }
        )

        if meta_abstract:

            resultat["Abstract"] = (
                meta_abstract.get(
                    "content",
                    ""
                )
            )


    return resultat


# ============================================================
# EXTRACTION GLOBALE
# ============================================================

def extraire_ers(nb_max):

    abstracts = []

    page = 0

    while len(abstracts) < nb_max:

        print()
        print(
            f"Page ERS {page + 1} "
            f"— {len(abstracts)}/{nb_max}"
        )


        resultats_page, page_suivante = (
            extraire_page(page)
        )


        if not resultats_page:

            print(
                "Aucun résultat trouvé."
            )

            break


        for item in resultats_page:

            if len(abstracts) >= nb_max:
                break


            print(
                f"[{len(abstracts) + 1}] "
                f"{item['Titre']}"
            )


            details = extraire_detail_abstract(
                item["URL"]
            )


            item.update(
                details
            )


            # =================================================
            # FILTRE SUR ANNÉE
            # =================================================

            date_pub = str(
                item.get(
                    "Date_publication",
                    ""
                )
            )


            conserver_annee = True

            if date_pub:

                try:

                    annee = int(
                        date_pub[:4]
                    )

                    conserver_annee = (
                        ANNEE_DEBUT
                        <= annee
                        <= ANNEE_FIN
                    )

                except:
                    pass


            if not conserver_annee:
                continue


            # =================================================
            # FILTRE RESPIRATOIRE
            # =================================================

            texte_complet = (
                str(item.get("Titre", ""))
                + " "
                + str(item.get("Abstract", ""))
                + " "
                + str(
                    item.get(
                        "Texte_resultat",
                        ""
                    )
                )
            )


            item["Respiratoire"] = (
                contient_un_terme(
                    texte_complet,
                    TERMES_RESPIRATOIRE
                )
            )


            # =================================================
            # FILTRE IA DE CONTRÔLE
            # =================================================

            item["IA"] = (
                contient_un_terme(
                    texte_complet,
                    TERMES_IA
                )
            )


            # On ne conserve que IA + respiratoire
            if (
                item["IA"]
                and
                item["Respiratoire"]
            ):

                abstracts.append(
                    item
                )


            time.sleep(
                PAUSE_REQUETES
            )


        if not page_suivante:

            print(
                "Dernière page atteinte."
            )

            break


        page += 1

        time.sleep(
            PAUSE_REQUETES
        )


    return abstracts


# ============================================================
# LANCEMENT
# ============================================================

abstracts = extraire_ers(
    NB_MAX_ABSTRACTS
)


print()
print(
    f"Nombre d'abstracts IA + respiratoire : "
    f"{len(abstracts)}"
)


# ============================================================
# DATAFRAME
# ============================================================

df = pd.DataFrame(
    abstracts
)


# ============================================================
# NETTOYAGE DES COLONNES
# ============================================================

if not df.empty:

    colonnes = [
        "Titre",
        "Auteurs",
        "Publication",
        "Date_publication",
        "DOI",
        "Abstract",
        "URL",
        "IA",
        "Respiratoire"
    ]


    colonnes_existantes = [
        colonne
        for colonne in colonnes
        if colonne in df.columns
    ]


    df = df[
        colonnes_existantes
    ]


    # --------------------------------------------------------
    # Suppression des doublons
    # --------------------------------------------------------

    if "DOI" in df.columns:

        df = df.drop_duplicates(
            subset=["DOI"],
            keep="first"
        )


    # --------------------------------------------------------
    # Tri par date
    # --------------------------------------------------------

    if "Date_publication" in df.columns:

        df["_date_tri"] = pd.to_datetime(
            df["Date_publication"],
            errors="coerce"
        )

        df = df.sort_values(
            "_date_tri",
            ascending=False
        )

        df = df.drop(
            columns=["_date_tri"]
        )


# ============================================================
# ONGLET PARAMÈTRES
# ============================================================

parametres = pd.DataFrame({

    "Parametre": [
        "ANNEE_DEBUT",
        "ANNEE_FIN",
        "NB_MAX_ABSTRACTS",
        "RESULTATS_PAR_PAGE",
        "TERMES_IA",
        "TERMES_RESPIRATOIRE"
    ],

    "Valeur": [
        ANNEE_DEBUT,
        ANNEE_FIN,
        NB_MAX_ABSTRACTS,
        RESULTATS_PAR_PAGE,
        " | ".join(TERMES_IA),
        " | ".join(
            TERMES_RESPIRATOIRE
        )
    ]
})


# ============================================================
# EXPORT EXCEL
# ============================================================

with pd.ExcelWriter(
    FICHIER_SORTIE,
    engine="openpyxl"
) as writer:

    df.to_excel(
        writer,
        sheet_name="ERS IA Respiratoire",
        index=False
    )

    parametres.to_excel(
        writer,
        sheet_name="Parametres",
        index=False
    )


print()
print("=" * 60)

print(
    "Extraction ERS terminée."
)

print(
    f"Fichier créé : "
    f"{FICHIER_SORTIE}"
)

print(
    f"Abstracts conservés : "
    f"{len(df)}"
)

print("=" * 60)