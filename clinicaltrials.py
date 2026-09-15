import requests
import pandas as pd
import time
import json

# ============================================================
# PARAMÈTRES
# ============================================================

# Nombre maximum d'essais à récupérer
NB_MAX_ESSAIS = 1000

# Taille des pages ClinicalTrials.gov
# Maximum autorisé par l'API : 1000
TAILLE_PAGE = 100

FICHIER_SORTIE = "C:/PYTHON/.data/clinicaltrials_IA_respiratoire.xlsx"

# ------------------------------------------------------------
# TERMES IA
# ------------------------------------------------------------

TERMES_IA = [
    '"artificial intelligence"',
    '"machine learning"',
    '"deep learning"',
    '"neural network"',
    '"neural networks"',
    '"large language model"',
    '"large language models"',
    '"foundation model"',
    '"computer vision"',
    '"natural language processing"',
    '"generative AI"',
    '"predictive model"'
]

# ------------------------------------------------------------
# TERMES RESPIRATOIRES
# ------------------------------------------------------------

TERMES_RESPIRATOIRE = [
    "respiratory",
    "pulmonary",
    "lung",
    "asthma",
    "COPD",
    '"chronic obstructive pulmonary disease"',
    "bronchiectasis",
    '"pulmonary fibrosis"',
    '"interstitial lung disease"',
    "spirometry",
    "emphysema",
    "airway",
    '"sleep apnea"',
    '"sleep apnoea"'
]


# ============================================================
# CONSTRUCTION DES REQUÊTES
# ============================================================

REQUETE_IA = " OR ".join(TERMES_IA)

REQUETE_RESPIRATOIRE = " OR ".join(TERMES_RESPIRATOIRE)

print("REQUÊTE IA :")
print(REQUETE_IA)

print()
print("REQUÊTE RESPIRATOIRE :")
print(REQUETE_RESPIRATOIRE)


# ============================================================
# SESSION HTTP
# ============================================================

session = requests.Session()

session.headers.update({
    "User-Agent": "ClinicalTrials research script"
})


# ============================================================
# FONCTIONS UTILES
# ============================================================

def get_nested(data, chemin, valeur_defaut=""):
    """
    Permet de récupérer facilement une donnée dans un JSON imbriqué.
    Exemple :
    get_nested(study, ["protocolSection", "identificationModule", "nctId"])
    """

    valeur = data

    for cle in chemin:

        if not isinstance(valeur, dict):
            return valeur_defaut

        valeur = valeur.get(cle)

        if valeur is None:
            return valeur_defaut

    return valeur


def liste_vers_texte(valeur):
    """
    Transforme une liste en texte lisible pour Excel.
    """

    if valeur is None:
        return ""

    if isinstance(valeur, list):

        valeurs = []

        for element in valeur:

            if isinstance(element, dict):
                valeurs.append(
                    json.dumps(
                        element,
                        ensure_ascii=False
                    )
                )

            else:
                valeurs.append(str(element))

        return " | ".join(valeurs)

    return str(valeur)


# ============================================================
# EXTRACTION DES INTERVENTIONS
# ============================================================

def extraire_interventions(study):

    interventions = get_nested(
        study,
        [
            "protocolSection",
            "armsInterventionsModule",
            "interventions"
        ],
        []
    )

    noms = []
    types = []
    descriptions = []

    for intervention in interventions:

        noms.append(
            intervention.get("name", "")
        )

        types.append(
            intervention.get("type", "")
        )

        descriptions.append(
            intervention.get("description", "")
        )

    return {
        "interventions": " | ".join(
            x for x in noms if x
        ),

        "types_interventions": " | ".join(
            x for x in types if x
        ),

        "descriptions_interventions": " | ".join(
            x for x in descriptions if x
        )
    }


# ============================================================
# EXTRACTION DES CRITÈRES DE JUGEMENT
# ============================================================

def extraire_outcomes(study):

    module = get_nested(
        study,
        [
            "protocolSection",
            "outcomesModule"
        ],
        {}
    )

    primaires = []
    secondaires = []

    for outcome in module.get("primaryOutcomes", []):

        mesure = outcome.get("measure", "")

        description = outcome.get(
            "description",
            ""
        )

        timeframe = outcome.get(
            "timeFrame",
            ""
        )

        texte = mesure

        if timeframe:
            texte += f" [{timeframe}]"

        if description:
            texte += f" : {description}"

        primaires.append(texte)


    for outcome in module.get(
        "secondaryOutcomes",
        []
    ):

        mesure = outcome.get("measure", "")

        description = outcome.get(
            "description",
            ""
        )

        timeframe = outcome.get(
            "timeFrame",
            ""
        )

        texte = mesure

        if timeframe:
            texte += f" [{timeframe}]"

        if description:
            texte += f" : {description}"

        secondaires.append(texte)


    return {
        "criteres_principaux": " | ".join(
            primaires
        ),

        "criteres_secondaires": " | ".join(
            secondaires
        )
    }


# ============================================================
# EXTRACTION DES LOCALISATIONS
# ============================================================

def extraire_localisations(study):

    lieux = get_nested(
        study,
        [
            "protocolSection",
            "contactsLocationsModule",
            "locations"
        ],
        []
    )

    pays = set()
    villes = set()
    centres = set()

    for lieu in lieux:

        pays.add(
            lieu.get("country", "")
        )

        villes.add(
            lieu.get("city", "")
        )

        centres.add(
            lieu.get("facility", "")
        )

    return {
        "pays": " | ".join(
            sorted(x for x in pays if x)
        ),

        "villes": " | ".join(
            sorted(x for x in villes if x)
        ),

        "centres": " | ".join(
            sorted(x for x in centres if x)
        )
    }


# ============================================================
# CONVERSION D'UNE ÉTUDE EN LIGNE EXCEL
# ============================================================

def convertir_etude(study):

    # --------------------------------------------------------
    # Identification
    # --------------------------------------------------------

    identification = get_nested(
        study,
        [
            "protocolSection",
            "identificationModule"
        ],
        {}
    )

    # --------------------------------------------------------
    # Statut
    # --------------------------------------------------------

    statut = get_nested(
        study,
        [
            "protocolSection",
            "statusModule"
        ],
        {}
    )

    # --------------------------------------------------------
    # Sponsor
    # --------------------------------------------------------

    sponsor_module = get_nested(
        study,
        [
            "protocolSection",
            "sponsorCollaboratorsModule"
        ],
        {}
    )

    sponsor = sponsor_module.get(
        "leadSponsor",
        {}
    )

    collaborateurs = sponsor_module.get(
        "collaborators",
        []
    )

    # --------------------------------------------------------
    # Design
    # --------------------------------------------------------

    design = get_nested(
        study,
        [
            "protocolSection",
            "designModule"
        ],
        {}
    )

    # --------------------------------------------------------
    # Conditions
    # --------------------------------------------------------

    conditions_module = get_nested(
        study,
        [
            "protocolSection",
            "conditionsModule"
        ],
        {}
    )

    # --------------------------------------------------------
    # Description
    # --------------------------------------------------------

    description = get_nested(
        study,
        [
            "protocolSection",
            "descriptionModule"
        ],
        {}
    )

    # --------------------------------------------------------
    # Éligibilité
    # --------------------------------------------------------

    eligibility = get_nested(
        study,
        [
            "protocolSection",
            "eligibilityModule"
        ],
        {}
    )


    # --------------------------------------------------------
    # Interventions
    # --------------------------------------------------------

    interventions = extraire_interventions(
        study
    )


    # --------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------

    outcomes = extraire_outcomes(
        study
    )


    # --------------------------------------------------------
    # Localisations
    # --------------------------------------------------------

    localisations = extraire_localisations(
        study
    )


    # --------------------------------------------------------
    # Collaborateurs
    # --------------------------------------------------------

    noms_collaborateurs = []

    for collaborateur in collaborateurs:

        nom = collaborateur.get(
            "name",
            ""
        )

        if nom:
            noms_collaborateurs.append(nom)


    # --------------------------------------------------------
    # Dates
    # --------------------------------------------------------

    start_date = statut.get(
        "startDateStruct",
        {}
    ).get("date", "")

    completion_date = statut.get(
        "completionDateStruct",
        {}
    ).get("date", "")

    primary_completion = statut.get(
        "primaryCompletionDateStruct",
        {}
    ).get("date", "")


    # --------------------------------------------------------
    # Nombre de participants
    # --------------------------------------------------------

    enrollment = design.get(
        "enrollmentInfo",
        {}
    )


    # --------------------------------------------------------
    # Création de la ligne
    # --------------------------------------------------------

    ligne = {

        "NCT_ID":
            identification.get(
                "nctId",
                ""
            ),

        "Titre":
            identification.get(
                "briefTitle",
                ""
            ),

        "Titre_officiel":
            identification.get(
                "officialTitle",
                ""
            ),

        "Statut":
            statut.get(
                "overallStatus",
                ""
            ),

        "Type_etude":
            design.get(
                "studyType",
                ""
            ),

        "Phases":
            liste_vers_texte(
                design.get(
                    "phases",
                    []
                )
            ),

        "Pathologies":
            liste_vers_texte(
                conditions_module.get(
                    "conditions",
                    []
                )
            ),

        "Mots_cles":
            liste_vers_texte(
                conditions_module.get(
                    "keywords",
                    []
                )
            ),

        "Resume":
            description.get(
                "briefSummary",
                ""
            ),

        "Description_detaillee":
            description.get(
                "detailedDescription",
                ""
            ),

        "Sponsor":
            sponsor.get(
                "name",
                ""
            ),

        "Type_sponsor":
            sponsor.get(
                "class",
                ""
            ),

        "Collaborateurs":
            " | ".join(
                noms_collaborateurs
            ),

        "Date_debut":
            start_date,

        "Date_fin_principale":
            primary_completion,

        "Date_fin":
            completion_date,

        "Nombre_participants":
            enrollment.get(
                "count",
                ""
            ),

        "Type_nombre_participants":
            enrollment.get(
                "type",
                ""
            ),

        "Sexe":
            eligibility.get(
                "sex",
                ""
            ),

        "Age_minimum":
            eligibility.get(
                "minimumAge",
                ""
            ),

        "Age_maximum":
            eligibility.get(
                "maximumAge",
                ""
            ),

        "Volontaires_sains":
            eligibility.get(
                "healthyVolunteers",
                ""
            ),

        "Criteres_inclusion_exclusion":
            eligibility.get(
                "eligibilityCriteria",
                ""
            ),

        "Interventions":
            interventions[
                "interventions"
            ],

        "Types_interventions":
            interventions[
                "types_interventions"
            ],

        "Descriptions_interventions":
            interventions[
                "descriptions_interventions"
            ],

        "Criteres_principaux":
            outcomes[
                "criteres_principaux"
            ],

        "Criteres_secondaires":
            outcomes[
                "criteres_secondaires"
            ],

        "Pays":
            localisations[
                "pays"
            ],

        "Villes":
            localisations[
                "villes"
            ],

        "Centres":
            localisations[
                "centres"
            ]
    }


    # --------------------------------------------------------
    # URL ClinicalTrials.gov
    # --------------------------------------------------------

    if ligne["NCT_ID"]:

        ligne["URL"] = (
            "https://clinicaltrials.gov/study/"
            + ligne["NCT_ID"]
        )

    else:

        ligne["URL"] = ""


    return ligne


# ============================================================
# EXTRACTION CLINICALTRIALS.GOV
# ============================================================

def extraire_essais(
    requete_ia,
    requete_respiratoire,
    nb_max_essais
):

    url = (
        "https://clinicaltrials.gov/api/v2/studies"
    )

    essais = []

    page_token = None

    while len(essais) < nb_max_essais:

        nombre_restant = (
            nb_max_essais - len(essais)
        )

        taille_page = min(
            TAILLE_PAGE,
            nombre_restant
        )

        params = {

            # Recherche IA dans les champs généraux
            "query.term":
                requete_ia,

            # Recherche respiratoire dans les champs
            # "Conditions or diseases"
            "query.cond":
                requete_respiratoire,

            "pageSize":
                taille_page,

            "format":
                "json"
        }


        if page_token:

            params["pageToken"] = (
                page_token
            )


        print()
        print(
            f"Extraction : "
            f"{len(essais)}/{nb_max_essais}"
        )


        try:

            response = session.get(
                url,
                params=params,
                timeout=90
            )

            response.raise_for_status()

        except requests.RequestException as e:

            print(
                f"Erreur API : {e}"
            )

            break


        data = response.json()

        studies = data.get(
            "studies",
            []
        )


        if not studies:

            print(
                "Aucun résultat supplémentaire."
            )

            break


        essais.extend(
            studies[:nombre_restant]
        )


        print(
            f"{len(studies)} essais récupérés"
        )


        # Jeton permettant de demander
        # la page suivante
        page_token = data.get(
            "nextPageToken"
        )


        if not page_token:

            print(
                "Dernière page atteinte."
            )

            break


        time.sleep(0.2)


    return essais[:nb_max_essais]


# ============================================================
# LANCEMENT DE L'EXTRACTION
# ============================================================

essais = extraire_essais(
    REQUETE_IA,
    REQUETE_RESPIRATOIRE,
    NB_MAX_ESSAIS
)

print()
print(
    f"Nombre total d'essais récupérés : "
    f"{len(essais)}"
)


# ============================================================
# CONVERSION EN DATAFRAME
# ============================================================

lignes = []

for compteur, study in enumerate(
    essais,
    start=1
):

    try:

        ligne = convertir_etude(
            study
        )

        lignes.append(
            ligne
        )

        print(
            f"[{compteur}/{len(essais)}] "
            f"{ligne['NCT_ID']} - "
            f"{ligne['Titre']}"
        )

    except Exception as e:

        print(
            f"Erreur sur l'essai "
            f"{compteur} : {e}"
        )


df = pd.DataFrame(
    lignes
)


# ============================================================
# TRI
# ============================================================

if not df.empty:

    if "Date_debut" in df.columns:

        df["_date_tri"] = pd.to_datetime(
            df["Date_debut"],
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
        "NB_MAX_ESSAIS",
        "TAILLE_PAGE",
        "REQUETE_IA",
        "REQUETE_RESPIRATOIRE"
    ],

    "Valeur": [
        NB_MAX_ESSAIS,
        TAILLE_PAGE,
        REQUETE_IA,
        REQUETE_RESPIRATOIRE
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
        sheet_name="Essais IA Respiratoire",
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
    "Extraction terminée."
)

print(
    f"Fichier créé : {FICHIER_SORTIE}"
)

print(
    f"Nombre d'essais : {len(df)}"
)

print("=" * 60)