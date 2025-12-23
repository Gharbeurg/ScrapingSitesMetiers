# Bibliotheque
import os
import time
import pandas as pd
import bs4
import re
import undetected_chromedriver as uc
from unidecode import unidecode

from selenium import webdriver
from datetime import datetime
from colorama import Fore, Back, Style

# variables
fichier_entree = "C:/DATA/github/.params/entree.txt"
url_pappers = "https://www.pappers.fr/recherche?q="
fichier_sortie_pappers = "C:/DATA/github/.data/pappers_content.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 10
temporisation = 3

df_posts = pd.DataFrame(columns =  ['societe', 'ape'])
df_posts = df_posts.reset_index(drop=True)
l_societe = []
l_ape = []

# suppression du fichier de sortie
print("{} - Suppression du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_sortie_pappers):
    os.remove(fichier_sortie_pappers)

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# lecture du fichier entree
print("{} - Lecture du fichier d'entree".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
with open(fichier_entree, 'r', encoding="utf-8", errors='ignore') as f:
    for societe in f:

        societe = unidecode(societe) #accents
        societe = societe.replace(' ', '+')
        societe = societe.strip()
  
        # Construction de la requête de recherche
        adresse = url_pappers + str(societe)
        driver.get(adresse)
        time.sleep(temporisation)
        
        #Recuperation du contenu
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')
        print(soup_contenu)

        # ape
        try:
            if soup_contenu.find("p", {"class":"key"}):
                liste_ape = soup_contenu.find_all("p", {"class":"key"})
                for ape in liste_ape:
                    code_ape = ape.find("p", {"class":"key"}).text
                    print(code_ape)
                    if 'NAF' in code_ape:
                        code_ape = code_ape[11:].strip()

            else:
                code_ape = ''

            #ajout dans la liste
            l_societe.append(societe)
            l_ape.append(code_ape)

            # parsing de chaque compte
            print("{} - société: {} - APE: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), societe, code_ape))

        except:
            print(Fore.RED + "Erreur sur cette societe." + Fore.RESET)

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['societe'] = l_societe
df_posts['ape'] = l_ape


# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_pappers):
    os.remove(fichier_sortie_pappers)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_pappers)

#fermeture
#driver.quit()