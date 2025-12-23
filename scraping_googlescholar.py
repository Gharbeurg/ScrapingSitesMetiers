#bibliothèques
import bs4
import time
import os
import pandas as pd
from datetime import datetime
from colorama import Fore, Back, Style

from selenium import webdriver
from selenium.webdriver.chrome.options import Options

#driver chrome
option = webdriver.ChromeOptions()
option.add_experimental_option('excludeSwitches', ['enable-logging'])
option.add_argument('--ignore-certificate-errors')
driver = webdriver.Chrome(options=option)

#variables
MAX_PAGE = 100
TERME_CHERCHE = 'apparel+industry' #AND OR
ANNEE_DEBUT='2021' #YYYY
temporisation = 5
nbre_erreur = 0

#requete scholar
query_scholar = "https://scholar.google.com/scholar?as_ylo={}&q={}&start={}"

fichier_sortie = "C:/DATA/github/.data/scholar_sortie.xlsx"
df_resultats = pd.DataFrame(columns =  ['titre', 'lien', 'auteurs', 'annee', 'journal','editeur', 'cite', 'resume'])
l_titre = []
l_lien = []
l_auteur = []
l_annee = []
l_journal = []
l_resume = []
l_cite = []
l_editeur =[]


# suppression du fichier de sortie
print("{} - Suppression du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_sortie):
    os.remove(fichier_sortie)

for numero_page in range(58, MAX_PAGE):
    try:
        # Numéro de page
        valeur_page = 10 * numero_page

        #construction de la requête
        query = query_scholar.format(ANNEE_DEBUT, TERME_CHERCHE,valeur_page)

        #recuperation du contenu de la page
        driver.get(query)
        get_url = driver.current_url
        time.sleep(temporisation)
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

        ##nombre de publications
        if numero_page == 0:
            nbre_publications = soup_contenu.find("div",  {"id":"gs_ab_md"}).text.strip()
            nbre_publications = nbre_publications.replace(" ", "")
            nbre_publications = nbre_publications.replace("Environ", "")
            position = nbre_publications.find("résultats")
            nbre_publications = nbre_publications[:position-1]

            print('{} - Nombre de publications : {}'.format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), nbre_publications))

        #Extraction du contenu de chaque résultat
        print("{} - Extraction des résultats de la page : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"),numero_page))
        result_block = soup_contenu.find_all('div', attrs = {'class': 'gs_ri'})
        for result in result_block:
            try:
                #extraction du lien
                result_url = result.find('a').get('href')

                #Extraction du titre
                if result.find("h3", {"class":"gs_rt"}):
                    result_titre = result.find("h3", {"class":"gs_rt"}).text
                else:
                    result_titre = "pas de titre"
            
                #Extraction du soustitre
                if result.find("div", {"class":"gs_a"}):
                    result_soustitre = result.find("div", {"class":"gs_a"}).text
                    position = result_soustitre.find("-")
                    result_auteurs = result_soustitre[:position-1]
                    result_soustitre = result_soustitre[position+1:]
                    position = result_soustitre.find("-")
                    result_editeur = result_soustitre[position+1:]
                    result_journal = result_soustitre[:position-1]
                    position = result_journal.find(",")
                    result_annee = result_journal[position+1:]
                    result_journal = result_journal[:position]
                else:
                    result_soustitre = "pas de soustitre"
                    result_auteurs = ""
                    result_journal = ""
                    result_annee = ""
                    result_editeur = ""

                #Extraction du contenu
                if result.find("div", {"class":"gs_rs"}):
                    result_description = result.find("div", {"class":"gs_rs"}).text
                else:
                    result_description = "pas de description"

                #Extraction des citations
                if result.find("div", {"class":"gs_fl"}):
                    result_citation = result.find("div", {"class":"gs_fl"}).text
                    result_citation = result_citation.replace("Enregistrer Citer", "")
                    result_citation = result_citation.replace(" Autres articles Les ", " - ")
                else:
                    result_citation = "pas de citation"
            
                print("{} - Etude : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"),result_url[0:50]))

            except:
                print(Fore.RED + "[+] Erreur sur cette publication : {}".format(result_url) + Fore.RESET)
                nbre_erreur +=1

            l_lien.append(result_url)
            l_titre.append(result_titre)
            l_auteur.append(result_auteurs)
            l_journal.append(result_journal)
            l_annee.append(result_annee)
            l_editeur.append(result_editeur)
            l_resume.append(result_description)
            l_cite.append(result_citation)

    except:
        print(Fore.RED + "[+] page introuvable" + Fore.RESET)

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_resultats['lien'] = pd.Series(l_lien)
df_resultats['titre'] = pd.Series(l_titre)
df_resultats['auteurs'] = pd.Series(l_auteur)
df_resultats['journal'] = pd.Series(l_journal)
df_resultats['annee'] = pd.Series(l_annee)
df_resultats['editeur'] = pd.Series(l_editeur)
df_resultats['resume'] = pd.Series(l_resume)
df_resultats['cite'] = pd.Series(l_cite)

# ecriture du fichier de résultats
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
with open(fichier_sortie, 'w', encoding="utf-8") as f:
    df_resultats.to_excel(fichier_sortie)
    
#fermeture du fichier
f.close()