#Bibliotheques
import pandas as pd
import bs4
import requests
import os
import re

from datetime import datetime
from colorama import Fore, Back, Style

#Variables
fichier_auteurs = "C:/DATA/github/.data/pubmed_auteurs.txt"
token_pubmed = "https://pubmed.ncbi.nlm.nih.gov/?term=" # adresse de pubmed
#filtre_annee = "&filter=years.2018-2021"
filtre_annee = ""
filtre_motscles = "clinical"
fichier_pubmed_sortie = "C:/DATA/github/.data/pubmed__nbre_article_sortie.xlsx"

df_articles = pd.DataFrame(columns =  ['auteurs','nb_articles'])
df_articles = df_articles.reset_index(drop=True)
l_auteurs = []
l_nb_articles = []
pages = []
nbre_pages_erreur = 0

# Collecte des pages à travailler
print("{} - Ouverture du fichier des auteurs".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
with open(fichier_auteurs, 'r', encoding="utf-8") as f:
    for line in f:
        try:
            # Remove all non-word characters (everything except numbers and letters)
            line = re.sub(r"[^\w\s]", '', line)

            # Replace all runs of whitespace with a plus
            line = re.sub(r"\s+", '+', line)

            #Minuscules
            line = line.lower()    

            # Construction de la requête
            requete = token_pubmed + filtre_motscles + line + filtre_annee
       
            # parsing de chaque auteur
            response = requests.get(requete)
            soup = bs4.BeautifulSoup(response.text, 'lxml')

            #Nombre articles, suppression des caractères inutiles
            nb_articles = re.sub(r"[^\w\s]", '', soup.find("div", {"class":"results-amount"}).text)
            nb_articles = re.sub('results', '', nb_articles)
            nb_articles = re.sub(',', '', nb_articles)
            nb_articles = re.sub(r"[^\w\s]", '', nb_articles)
            nb_articles = re.sub(r"\s+", '', nb_articles)
            if nb_articles == 'Nowerefound':
                nb_articles = 0
            
            l_nb_articles.append(nb_articles)

            #Auteur
            l_auteurs.append(line)
            print("{} - Traitement auteur : {} - Articles : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), line, nb_articles))
        
        except:
            print("{} - Traitement auteur : {} - Articles : 1".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), line))
            l_nb_articles.append('1')
            l_auteurs.append(line)

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles['nb_articles'] = l_nb_articles
df_articles['auteurs'] = l_auteurs

# ecriture des fichiers de sortie
if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles.to_excel(fichier_pubmed_sortie)

# Fermture du fichier des auteurs
print("{} - Fermeture du fichier des auteurs".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f.close()
