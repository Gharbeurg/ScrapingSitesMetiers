#bibliothèques
import os
import pandas as pd
import re
import time
import Bio


from Bio import Entrez
from Bio import Medline
from datetime import datetime
from colorama import Fore, Back, Style

#variables
MAX_COUNT = 99
TERME_CHERCHE = 'leukemia artificial intelligence' #AND OR
DATE_DEBUT='2025/01/01' #YYYY/M
DATE_FIN='2025/11/01' #YYYY/M
MAIL_PUBMED = 'dev@atawao.com'
BASE_NCBI = 'pubmed'
ORDRE_TRI = 'pub date' #pub date, relevance, first author, last author, title, journal
FORMAT_RESULTAT = 'medline' #XML, medline, uilist pour les id, abstract
FORMAT_ENREGISTREMENT = 'text' #XML ou Texte

#API KEY CONFIGURATION
APIKEY = "2bad5e3552250d73f6000974773fb3d05c08"
os.environ['NCBI_API_KEY'] = APIKEY

#bibliotheque a importer après la définition de la variable dans l'environnement
import metapub

fichier_facteur_impact = "C:/DATA/code/.params/impact_factor.txt"
fichier_pubmed_sortie = "C:/DATA/code/.data/pubmed_sortie.xlsx"
df_pages = pd.DataFrame(columns =  ['id', 'titre', 'Ifactor', 'DOI', 'PMC', 'OID', 'OCI', 'copyright', 'auteurs', 'journal', 'creation', 'resume'])
l_id = []
l_titre = []
l_factor = []
l_doi = []
l_pmc = []
l_oid = []
l_oci = []
l_copyright = []
l_auteurs = []
l_journal = []
l_creation = []
l_resume =[]

#Ouverture du fichier d'entrée
print("{} - Ouverture du fichier des facteurs d impact".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_impact = pd.read_csv(fichier_facteur_impact, sep=';', header=0)
df_impact.columns = ["index", "titre", "if"]

#construction de la requête
TERM = f"{TERME_CHERCHE} AND ({DATE_DEBUT}[PDAT] : {DATE_FIN}[PDAT])"

Entrez.email = MAIL_PUBMED
Entrez.api_key = os.environ['NCBI_API_KEY']
h = Entrez.esearch(db=BASE_NCBI, retmax=MAX_COUNT, term=TERM) #requete
result = Entrez.read(h)
ids = result['IdList']

print('{} - Nombre de publications avec le terme {}: {}'.format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), TERME_CHERCHE, result['Count']))

#Lecture de chaque enregistrement et récupération du détail de l'article
h = Entrez.efetch(db=BASE_NCBI, id=ids, rettype=FORMAT_RESULTAT, retmode=FORMAT_ENREGISTREMENT)
records = Medline.parse(h)

#affichage des résultats
for record in records:

    try:
        time.sleep(1)  # Attendre 1 seconde entre chaque requête
        #recherche du facteur impact
        journal = record.get('JT', '?')
        journal = journal.lower()

        # Check if a value exists in the 'Name' column
        if journal in df_impact['titre'].values:
            facteur_impact=df_impact.loc[df_impact['titre'] == journal, 'if'].iloc[0]
        else:
            facteur_impact = ''

        if metapub.FindIt(record.get('PMID', '?')).doi:
            lien_DOI = metapub.FindIt(record.get('PMID', '?')).doi
        else :
            lien_DOI = ""
        
        l_id.append(record.get('PMID', '?'))
        l_titre.append(record.get('TI', '?'))
        l_factor.append(facteur_impact)
        l_doi.append(lien_DOI)
        l_copyright.append(record.get('CI', '?'))
        l_oid.append(record.get('OID', '?'))
        l_pmc.append(record.get('PMC', '?'))
        l_oci.append(record.get('OCI', '?'))
        l_auteurs.append(record.get('AU', '?'))
        l_journal.append(record.get('JT', '?'))
        l_creation.append(record.get('DEP', '?'))
        l_resume.append(record.get('AB', '?'))

        #autres champs disponibles
        #date_derniere_revision = record.get('LR', '?')
        #date_creation = record.get('CRDT', '?')
        #titre_journal = record.get('JT', '?')
        #auteur_identifiant = record.get('AUID', '?')
        #langue = record.get('LA', '?')
        #autre_resume = record.get('OAB', '?')
        #statut_publication = record.get('PST', '?')
        #source = record.get('SO', '?')

    except:
        print(Fore.RED + "{} - ca déconne : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), record.get('DEP', '?')) + Fore.RESET)

#creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_pages['id'] = l_id
df_pages['titre'] = pd.Series(l_titre)
df_pages['Ifactor'] = pd.Series(l_factor)
df_pages['DOI'] = pd.Series(l_doi)
df_pages['PMC'] = pd.Series(l_pmc)
df_pages['OID'] = pd.Series(l_oid)
df_pages['OCI'] = pd.Series(l_oci)
df_pages['copyright'] = pd.Series(l_copyright)
df_pages['auteurs'] = pd.Series(l_auteurs)
df_pages['journal'] = pd.Series(l_journal)
df_pages['creation'] = pd.Series(l_creation)
df_pages['resume'] = pd.Series(l_resume)

# ecriture des fichiers de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

df_pages.to_excel(fichier_pubmed_sortie)