#bibliotheques
import pandas as pd
import unicodedata
import bs4
import requests
import os
from tabulate import tabulate
from datetime import datetime

#variables
fichier_sujet = "C:/DATA/github/.data/forum_sujets.xlsx"
token_forum = "https://forum.doctissimo.fr/sante/arthrose-os/liste_sujet-1.htm" # adresse de la premiere page du forum
nbre_pages_forum = 20
nbre_pages_erreur = 0
df_sujets = pd.DataFrame(columns =  ['sujet','auteur','Rep','lus'])
df_sujets = df_sujets.reset_index(drop=True)
l_sujet = []
l_auteur = []
l_reponse = []
l_lu = []

# collecte des pages du forum
def get_pages(token, nb):
    pages = []
    for i in range(1,nb+1):
        j = token + str(i) + '.htm'
        pages.append(j)
    return pages

print("{} - Collecte des pages du forum".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
pages = get_pages(token_forum,nbre_pages_forum)

# parsing de chaque page sujet
for i in pages:
    try:
        response = requests.get(i)
        #time.sleep(random.randrange(1,5))
        print("[+] traitement de la page sujet : ", i)
        soup = bs4.BeautifulSoup(response.text, 'lxml')

        tag = soup.find_all("td", {"class":"sujetCase3"})
        for el in tag:
            chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
            chaine = str(chaine, "utf-8").lower()
            l_sujet.append(chaine)

        tag = soup.find_all("td", {"class":"sujetCase6"})
        for el in tag:
            chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
            chaine = str(chaine, "utf-8").lower()
            l_auteur.append(chaine)

        tag = soup.find_all("td", {"class":"sujetCase7"})
        for el in tag:
           chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
           chaine = str(chaine, "utf-8").lower().lower().replace(" ","")
           chaine = int(chaine)
           l_reponse.append(chaine)

        tag = soup.find_all("td", {"class":"sujetCase8"})
        for el in tag:
           chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
           chaine = str(chaine, "utf-8").lower().replace(" ","")
           chaine = int(chaine)
           l_lu.append(chaine)

    except:
        nbre_pages_erreur += 1

print("[+] Nombe de pages sujets en erreur : ", nbre_pages_erreur )
nbre_pages_erreur = 0

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_sujets['sujet'] = l_sujet
df_sujets['auteur'] = l_auteur
df_sujets['Rep'] = l_reponse
df_sujets['lus'] = l_lu

# ecriture du fichier de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_sujet):
    os.remove(fichier_sujet)
df_sujets.to_excel(fichier_sujet)