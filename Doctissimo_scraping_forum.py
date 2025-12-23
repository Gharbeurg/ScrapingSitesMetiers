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
fichier_adresses = "C:/DATA/github/.data/forum_adresse.xls"
fichier_contenu = "C:/DATA/github/.data/forum_contenu.txt"
token_forum = "https://forum.doctissimo.fr/sante/cholesterol/liste_sujet-"
nbre_pages_forum = 7
nbre_pages_erreur = 0
df_sujets = pd.DataFrame(columns =  ['sujet','auteur','Rep','lus'])
df_sujets = df_sujets.reset_index(drop=True)
df_adresses = pd.DataFrame(columns =  ['lien', 'sujet'])
l_lien = []
l_sujet = []
l_auteur = []
l_reponse = []
l_lu = []
l_sujetlien = []
l_contenu = []

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
        print("[+] Traitement de la page sujet : ", i)
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

# parsing de chaque page lien vers le detail
for i in pages:
    try:
        response = requests.get(i)
        print("[+] collecte des liens vers la page de détail :", i)
        soup = bs4.BeautifulSoup(response.text, 'lxml')

        tag = soup.find_all("a", {"class":"cCatTopic"})
        for el in tag:
            l_lien.append(el.get('href'))
            l_sujetlien.append(el.text)

    except:
        nbre_pages_erreur += 1

print("[+] Nombe de pages adresses en erreur : ", nbre_pages_erreur )
nbre_pages_erreur = 0

# parsing de chaque page de détail
for i in l_lien:
    try:
        response = requests.get(i)
        print("[+] parsing de la page de détail :", i)
        soup = bs4.BeautifulSoup(response.text, 'lxml')

        tag = soup.find_all("span", {"itemprop": "text", "hidden": True})
        for el in tag:
            chaine = unicodedata.normalize('NFKD', el.text).encode('ASCII', 'ignore')
            chaine = str(chaine, "utf-8").lower()
            l_contenu.append(chaine)

    except:
        print("[+] erreur pour cette page")

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_sujets['sujet'] = l_sujet
df_sujets['auteur'] = l_auteur
df_sujets['Rep'] = l_reponse
df_sujets['lus'] = l_lu
df_adresses['lien'] = l_lien
df_adresses['sujet'] = l_sujetlien

# ecriture du fichier sujet
print("{} - Ecriture des fichiers de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_sujet):
    os.remove(fichier_sujet)
df_sujets.to_excel(fichier_sujet)

# ecriture du fichier adresse
if os.path.exists(fichier_adresses):
    os.remove(fichier_adresses)
df_adresses.to_excel(fichier_adresses)

# ecriture du fichier contenu
if os.path.exists(fichier_contenu):
    os.remove(fichier_contenu)
with open(fichier_contenu, 'w') as f:
    for item in l_contenu:
        # item = item.encode('utf8')
        f.write(f'{item}\n')