# Bibliotheque
import os
import time
import pandas as pd
import bs4
import re

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from nordvpn_switcher import initialize_VPN, rotate_VPN
from datetime import datetime
from colorama import Fore, Back, Style

# variables
url_bpi = "https://lehub.web.bpifrance.fr/search?refinementList[thematique][0]=Health%20Tech%20/%20Medtech&page="
fichier_sortie_bpi = "C:/DATA/github/.data/bpi_content.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 105
temporisation = 4

df_posts = pd.DataFrame(columns =  ['startup_lien', 'startup_nom'])
df_posts = df_posts.reset_index(drop=True)
l_startup_lien = []
l_startup_nom = []

# OUverture d'un VPN
print("{} - Ouverture du VPN".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
settings = initialize_VPN(area_input=['Belgium,France,Netherlands,Germany,Slovenia,Switzerland,Denmark,Croatia,Finland,Estonia,Portugal,Greece, Serbia,Austria,Spain,Poland'])
rotate_VPN(settings) #refer to the instructions variable here
time.sleep(temporisation)

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# ouverture du driver
driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# Collecte des liens
print("{} - Démarrage de la collecte des liens".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

for i in range(1, page_max):
    
    # parsing de chaque compte
    print("{} - page {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i))

    # Construction de la requête de recherche
    adresse = url_bpi + str(i)
    driver.get(adresse)
    time.sleep(temporisation)

    #Recuperation du contenu
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    # liste des liens
    liste_startups = soup_contenu.find_all("a", {"class":"sc-jKTccl cnsQRA"})
    for startup in liste_startups:
        try:
            # lien ver la fiche startup
            lien_fiche_startup = startup.get('href')
            lien_fiche_startup = 'https://lehub.web.bpifrance.fr' + lien_fiche_startup

            #ajout dans la liste
            l_startup_lien.append(lien_fiche_startup)
            l_startup_nom.append('')

        except:
            print(Fore.RED + "Erreur sur cette page." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['startup_lien'] = l_startup_lien
df_posts['startup_nom'] = l_startup_nom

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_bpi):
    os.remove(fichier_sortie_bpi)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_bpi)

#fermeture
driver.quit()