#bibliotheques
import pandas as pd
import time
import bs4
import random
import requests
import os
import pickle
from tabulate import tabulate

#variables
fichier_annuaire = "C:/DATA/github/.data/splf_annuaire.xls"
token_forum = "http://ig.splf.fr/?p="
nbre_pages_forum = 30
nbre_pages_erreur = 0
df_annuaire = pd.DataFrame(columns =  ['nom','adresse'])
df_annuaire = df_annuaire.reset_index(drop=True)
l_nom = []
l_adresse = []

# collecte des pages du forum
def get_pages(token, nb):
    pages = []
    for i in range(1,nb+1):
        j = token + str(i) + '.htm'
        pages.append(j)
    return pages

pages = get_pages(token_forum,nbre_pages_forum)

# parsing de chaque page sujet
for i in pages:
    try:
        response = requests.get(i)
        #time.sleep(random.randrange(1,5))
        soup = bs4.BeautifulSoup(response.text, 'lxml')

        tag = tag = soup.find_all("td", {"width": "45%"})
        for el in tag:
            chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
            chaine = str(chaine, "utf-8").lower()
            l_nom.append(chaine)

        tag = soup.find_all("td", {"width": "45%"})
        for el in tag:
            chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
            chaine = str(chaine, "utf-8").lower()
            l_adresse.append(chaine)

    except:
        nbre_pages_erreur += 1

print("Nombe de pages sujets en erreur : ", nbre_pages_erreur )
nbre_pages_erreur = 0

# creation du dataframe
df_annuaire['nom'] = l_nom
#df_annuaire['adresse'] = l_adresse

# ecriture du fichier sujet
if os.path.exists(fichier_annuaire):
    os.remove(fichier_annuaire)
df_annuaire.to_excel(fichier_annuaire)