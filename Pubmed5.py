#bibliothèques
from Bio import Entrez
from Bio import Medline
import os
import pandas as pd
import re
import bs4

from datetime import datetime

#variables
MAX_COUNT = 3000
TERME_CHERCHE = 'chatgpt' #AND OR
DATE_DEBUT='2023/1/1' #YYYY/M
DATE_FIN='2024/2/1' #YYYY/M
MAIL_PUBMED = 'dev@atawao.com'
BASE_NCBI = 'pubmed'
ORDRE_TRI = 'pub date' #pub date, relevance, first author, last author, title, journal
FORMAT_RESULTAT = 'medline' #XML, medline, uilist pour les id, abstract
FORMAT_ENREGISTREMENT = 'text' #XML ou Texte

fichier_facteur_impact = "C:/PYTHON/.params/impact_factor.txt"
fichier_pubmed_sortie = "C:/PYTHON/.data/pubmed_sortie.xlsx"
df_pages = pd.DataFrame(columns =  ['id', 'titre', 'Ifactor', 'auteurs', 'journal', 'creation', 'resume'])
l_id = []
l_titre = []
l_factor = []
l_auteurs = []
l_journal = []
l_creation = []
l_resume =[]

#Ouverture du fichier d'entrée
print("{} - Ouverture du fichier des facteurs d impact".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_impact = pd.read_csv(fichier_facteur_impact, sep=';', header=0)
df_impact.columns = ["index", "titre", "if"]

#construction de la requête
TERM = TERME_CHERCHE + ' ' + DATE_DEBUT + ':' + DATE_FIN + '[Publication Date]'

Entrez.email = MAIL_PUBMED
h = Entrez.esearch(db=BASE_NCBI, retmax=MAX_COUNT, term=TERM) #requete
result = Entrez.read(h)
ids = result['IdList']

print('{} - Nombre de publications avec le terme {}: {}'.format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), TERME_CHERCHE, result['Count']))

#Lecture de chaque enregitrement et récupération du détail de l'article
h = Entrez.efetch(db=BASE_NCBI, id=ids, rettype=FORMAT_RESULTAT, retmode=FORMAT_ENREGISTREMENT)
records = Medline.parse(h)

#affichage des résultats
for record in records:

    #recherche du facteur impact
    journal = record.get('JT', '?')
    journal = journal.lower()

    # Check if a value exists in the 'Name' column
    if journal in df_impact['titre'].values:
        facteur_impact=df_impact.loc[df_impact['titre'] == journal, 'if'].iloc[0]
    else:
        facteur_impact = ''

    l_id.append(record.get('PMID', '?'))
    l_titre.append(record.get('TI', '?'))
    l_factor.append(facteur_impact)
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

#creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_pages['id'] = l_id
df_pages['titre'] = pd.Series(l_titre)
df_pages['Ifactor'] = pd.Series(l_factor)
df_pages['auteurs'] = pd.Series(l_auteurs)
df_pages['journal'] = pd.Series(l_journal)
df_pages['creation'] = pd.Series(l_creation)
df_pages['resume'] = pd.Series(l_resume)

# ecriture des fichiers de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

df_pages.to_excel(fichier_pubmed_sortie)