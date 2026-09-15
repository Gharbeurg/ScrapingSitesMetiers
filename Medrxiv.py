import requests
import pandas as pd
import time
import xml.etree.ElementTree as ET
from urllib.parse import urljoin

# ============================================================
# PARAMÈTRES
# ============================================================

DATE_DEBUT = "2026-01-01"
DATE_FIN = "2026-09-10"

# Nombre maximum d'articles medRxiv à récupérer AVANT filtrage
NB_MAX_ARTICLES = 1000

FICHIER_SORTIE = "C:/PYTHON/.data/medrxiv_IA_respiratoire_2026.xlsx"

# Pause entre deux téléchargements de texte intégral
PAUSE_XML = 0.2

TERMES_IA = [
    "artificial intelligence",
    "machine learning",
    "deep learning",
    "neural network",
    "neural networks",
    "large language model",
    "large language models",
    "llm",
    "foundation model",
    "transformer",
    "computer vision"
]

TERMES_RESPIRATOIRE = [
    "respiratory",
    "pulmonary",
    "lung",
    "lungs",
    "asthma",
    "copd",
    "chronic obstructive pulmonary disease",
    "bronchial",
    "bronchitis",
    "bronchiectasis",
    "pneumonia",
    "pulmonary fibrosis",
    "interstitial lung disease",
    "spirometry",
    "spirometric",
    "emphysema",
    "pleural",
    "airway",
    "airways",
    "respiration"
]


# ============================================================
# SESSION HTTP
# ============================================================

session = requests.Session()

session.headers.update({
    "User-Agent": "Mozilla/5.0 medRxiv research script"
})


# ============================================================
# 1. EXTRACTION MEDRXIV AVEC PAGINATION
# ============================================================

def extraire_medrxiv(date_debut, date_fin, nb_max_articles):

    articles = []
    cursor = 0

    while len(articles) < nb_max_articles:

        url = (
            f"https://api.medrxiv.org/details/medrxiv/"
            f"{date_debut}/{date_fin}/{cursor}/json"
        )

        print(
            f"Extraction medRxiv : {len(articles)}/{nb_max_articles}"
        )

        try:
            response = session.get(url, timeout=60)
            response.raise_for_status()

        except requests.RequestException as e:
            print(f"Erreur API medRxiv : {e}")
            break

        data = response.json()

        collection_originale = data.get("collection", [])

        if not collection_originale:
            break

        nb_restant = nb_max_articles - len(articles)

        collection = collection_originale[:nb_restant]

        articles.extend(collection)

        print(
            f"  {len(collection)} articles ajoutés "
            f"— total : {len(articles)}"
        )

        if len(articles) >= nb_max_articles:
            break

        # IMPORTANT :
        # on teste la taille de la réponse originale de l'API
        if len(collection_originale) < 100:
            break

        cursor += 100

        time.sleep(0.2)

    return articles


# ============================================================
# 2. FILTRAGE PAR TERMES
# ============================================================

def contient_terme(texte, liste_termes):

    texte = str(texte).lower()

    return any(
        terme.lower() in texte
        for terme in liste_termes
    )


# ============================================================
# 3. NETTOYAGE DU TEXTE XML
# ============================================================

def texte_element(element):

    if element is None:
        return ""

    morceaux = []

    for texte in element.itertext():

        texte = texte.strip()

        if texte:
            morceaux.append(texte)

    return " ".join(morceaux)


# ============================================================
# 4. IDENTIFICATION DES SECTIONS
# ============================================================

def classer_section(titre):

    titre = titre.lower().strip()

    # Méthodes
    if any(x in titre for x in [
        "method",
        "methods",
        "methodology",
        "materials and methods",
        "study design",
        "experimental procedures"
    ]):
        return "methodes"

    # Résultats
    if any(x in titre for x in [
        "result",
        "results",
        "findings"
    ]):
        return "resultats"

    # Discussion
    if any(x in titre for x in [
        "discussion"
    ]):
        return "discussion"

    # Conclusion
    if any(x in titre for x in [
        "conclusion",
        "conclusions"
    ]):
        return "conclusion"

    # Introduction
    if any(x in titre for x in [
        "introduction",
        "background"
    ]):
        return "introduction"

    return None


# ============================================================
# 5. EXTRACTION DES SECTIONS DU XML
# ============================================================

def extraire_sections_xml(contenu_xml):

    resultat = {
        "texte_complet": "",
        "introduction": "",
        "methodes": "",
        "resultats": "",
        "discussion": "",
        "conclusion": ""
    }

    try:
        root = ET.fromstring(contenu_xml)

    except ET.ParseError:
        return resultat

    # --------------------------------------------------------
    # Texte complet du corps de l'article
    # --------------------------------------------------------

    body = root.find(".//body")

    if body is not None:
        resultat["texte_complet"] = texte_element(body)

    else:
        resultat["texte_complet"] = texte_element(root)

    # --------------------------------------------------------
    # Sections structurées
    # --------------------------------------------------------

    for section in root.findall(".//sec"):

        titre_element = section.find("./title")

        if titre_element is None:
            continue

        titre = texte_element(titre_element)

        categorie = classer_section(titre)

        if categorie is None:
            continue

        contenu = texte_element(section)

        # Évite les doublons issus de sous-sections
        if contenu and contenu not in resultat[categorie]:

            if resultat[categorie]:
                resultat[categorie] += "\n\n"

            resultat[categorie] += contenu

    return resultat


# ============================================================
# 6. CONSTRUCTION DE L'URL XML
# ============================================================

def construire_url_xml(xml_path):

    if pd.isna(xml_path):
        return None

    xml_path = str(xml_path).strip()

    if not xml_path:
        return None

    # Cas où l'API fournit directement une URL complète
    if xml_path.startswith("http://") or xml_path.startswith("https://"):
        return xml_path

    # Cas où l'API fournit un chemin relatif
    return urljoin(
        "https://www.medrxiv.org/",
        xml_path
    )


# ============================================================
# 7. TÉLÉCHARGEMENT DU TEXTE COMPLET
# ============================================================

def recuperer_texte_integral(xml_path):

    url_xml = construire_url_xml(xml_path)

    resultat = {
        "url_xml": url_xml,
        "statut_xml": "",
        "texte_complet": "",
        "introduction": "",
        "methodes": "",
        "resultats": "",
        "discussion": "",
        "conclusion": ""
    }

    if not url_xml:

        resultat["statut_xml"] = "XML non indiqué par l'API"

        return resultat

    try:

        response = session.get(
            url_xml,
            timeout=60
        )

        response.raise_for_status()

        sections = extraire_sections_xml(
            response.content
        )

        resultat.update(sections)

        if resultat["texte_complet"]:

            resultat["statut_xml"] = "OK"

        else:

            resultat["statut_xml"] = "XML téléchargé mais texte non extrait"

    except requests.RequestException as e:

        resultat["statut_xml"] = (
            f"Erreur téléchargement : {str(e)}"
        )

    except Exception as e:

        resultat["statut_xml"] = (
            f"Erreur traitement XML : {str(e)}"
        )

    return resultat


# ============================================================
# 8. EXTRACTION DES MÉTADONNÉES
# ============================================================

articles = extraire_medrxiv(
    DATE_DEBUT,
    DATE_FIN,
    NB_MAX_ARTICLES
)

df = pd.DataFrame(articles)

print()
print(
    f"Nombre total d'articles medRxiv récupérés : {len(df)}"
)


# ============================================================
# 9. FILTRAGE IA + RESPIRATOIRE
# ============================================================

if not df.empty:

    df["texte_analyse"] = (
        df["title"].fillna("").astype(str)
        + " "
        + df["abstract"].fillna("").astype(str)
    ).str.lower()

    df["IA"] = df["texte_analyse"].apply(
        lambda x: contient_terme(
            x,
            TERMES_IA
        )
    )

    df["Respiratoire"] = df["texte_analyse"].apply(
        lambda x: contient_terme(
            x,
            TERMES_RESPIRATOIRE
        )
    )

    resultats = df[
        df["IA"] & df["Respiratoire"]
    ].copy()

    print(
        f"Articles IA + respiratoire identifiés : "
        f"{len(resultats)}"
    )


    # ========================================================
    # 10. RÉCUPÉRATION DU TEXTE INTÉGRAL
    # ========================================================

    colonnes_xml = [
        "url_xml",
        "statut_xml",
        "texte_complet",
        "introduction",
        "methodes",
        "resultats",
        "discussion",
        "conclusion"
    ]

    for colonne in colonnes_xml:
        resultats[colonne] = ""


    for compteur, index in enumerate(
        resultats.index,
        start=1
    ):

        titre = resultats.at[index, "title"]

        print()
        print(
            f"[{compteur}/{len(resultats)}] "
            f"{titre}"
        )

        xml_path = resultats.at[
            index,
            "jats xml path"
        ] if "jats xml path" in resultats.columns else None

        donnees_xml = recuperer_texte_integral(
            xml_path
        )

        for colonne, valeur in donnees_xml.items():

            resultats.at[
                index,
                colonne
            ] = valeur

        print(
            f"  XML : "
            f"{donnees_xml['statut_xml']}"
        )

        time.sleep(PAUSE_XML)


    # ========================================================
    # 11. TRI
    # ========================================================

    if "date" in resultats.columns:

        resultats["date"] = pd.to_datetime(
            resultats["date"],
            errors="coerce"
        )

        resultats = resultats.sort_values(
            by="date",
            ascending=False
        )


    # ========================================================
    # 12. COLONNES DE L'ONGLET PRINCIPAL
    # ========================================================

    colonnes_souhaitees = [
        "date",
        "title",
        "authors",
        "abstract",
        "category",
        "doi",
        "version",
        "type",
        "license",
        "published",

        "statut_xml",
        "url_xml",

        "introduction",
        "methodes",
        "resultats",
        "discussion",
        "conclusion",

        "texte_complet"
    ]

    colonnes_existantes = [
        colonne
        for colonne in colonnes_souhaitees
        if colonne in resultats.columns
    ]

    resultats_export = resultats[
        colonnes_existantes
    ].copy()


    # ========================================================
    # 13. STATISTIQUES
    # ========================================================

    nb_xml_ok = (
        resultats["statut_xml"] == "OK"
    ).sum()

    print()
    print(
        f"Textes intégraux récupérés : "
        f"{nb_xml_ok}/{len(resultats)}"
    )


    # ========================================================
    # 14. EXPORT EXCEL
    # ========================================================

    with pd.ExcelWriter(
        FICHIER_SORTIE,
        engine="openpyxl"
    ) as writer:

        # Articles pertinents + texte intégral
        resultats_export.to_excel(
            writer,
            sheet_name="IA Respiratoire",
            index=False
        )

        # Tous les articles avant filtrage
        df.drop(
            columns=["texte_analyse"],
            errors="ignore"
        ).to_excel(
            writer,
            sheet_name="Tous les articles",
            index=False
        )


    print()
    print("Extraction terminée.")
    print(
        f"Fichier créé : {FICHIER_SORTIE}"
    )


else:

    print("Aucun article récupéré.")