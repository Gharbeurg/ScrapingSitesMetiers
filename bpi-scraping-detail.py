# Bibliotheque
import os
import time
import pandas as pd
import bs4
import re

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from nordvpn_switcher import initialize_VPN, rotate_VPN
from datetime import datetime
from colorama import Fore, Back, Style

# variables
fichier_entree_bpi = "C:/DATA/github/.data/bpi_startups.txt"
fichier_sortie_bpi = "C:/DATA/github/.data/bpi_content_detail.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 105
temporisation = 4

df_posts = pd.DataFrame(columns =  ['startup_nom', 'startup_creation', 'startup_montant', 'startup_effectif', 'startup_identite', 'startup_contact', 'startup_secteur', 'startup_metier', 'startup_objectif', 'startup_thematique', 'startup_typeproduit', 'startup_business_model', 'startup_accelerateur', 'startup_presence', 'startup_clients', 'startup_produits', 'startup_lien'])
df_posts = df_posts.reset_index(drop=True)
l_startup_lien = []
l_startup_nom = []
l_startup_creation = []
l_startup_montant = []
l_startup_effectif = []
l_startup_identite = []
l_startup_contact = []
l_startup_secteur  = []
l_startup_metier  = []
l_startup_objectif  = []
l_startup_thematique  = []
l_startup_typeproduit  = []
l_startup_business_model  = []
l_startup_accelerateur  = []
l_startup_presence  = []
l_startup_clients = []
l_startup_produits = []

# OUverture d'un VPN
print("{} - Ouverture du VPN".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
settings = initialize_VPN(area_input=['Belgium,France,Netherlands,Germany,Slovenia,Switzerland,Denmark,Croatia,Finland,Estonia,Portugal,Greece, Serbia,Austria,Spain,Poland'])
rotate_VPN(settings) #refer to the instructions variable here
time.sleep(temporisation)

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
option = webdriver.ChromeOptions()
driver = webdriver.Chrome(options = option)
#driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# parsing de chaque startup
print("{} - Parsing de chaque startup".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

with open(fichier_entree_bpi, 'r', encoding="utf-8") as f:
    for line in f:
        nombre_requete +=1

        driver.get(line)
        time.sleep(temporisation)

        #Recuperation du contenu
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

        try:
            # Nom de la societe
            if soup_contenu.find('h1', attrs = {'class': 'sc-fTQvRK dgiayt'}):
                startup_nom = soup_contenu.find('h1', attrs = {'class': 'sc-fTQvRK dgiayt'}).text
            else : startup_nom = ''
    
            # Description de la societe
            if soup_contenu.find('p', attrs = {'class': 'sc-bUhFKy bgkZQl'}):
                startup_description = soup_contenu.find('p', attrs = {'class': 'sc-bUhFKy bgkZQl'}).text
            else : startup_description = ''
      
            # creation, montant, effectifs de la societe
            if soup_contenu.find('div', attrs = {'class': 'sc-bjeSbO iIFufF'}):
                startup_creation = soup_contenu.find('div', attrs = {'class': 'sc-bjeSbO iIFufF'}).text
            else : startup_creation = ''
        
            # identité de la societe
            if soup_contenu.find('div', attrs = {'class': 'sc-dkqQuH knuDIX startup__identity'}):
                startup_identite = soup_contenu.find('div', attrs = {'class': 'sc-dkqQuH knuDIX startup__identity'}).text

            # clients de la societe
            if soup_contenu.find('div', attrs = {'class': 'sc-dwsnSq ecpChw startup__clients'}):
                startup_clients = soup_contenu.find('div', attrs = {'class': 'sc-dwsnSq ecpChw startup__clients'}).text
            else : startup_clients = ''

            # produits de la societe
            if soup_contenu.find('div', attrs = {'class': 'sc-dwsnSq ecpChw startup__products'}):
                startup_produits = soup_contenu.find('div', attrs = {'class': 'sc-dwsnSq ecpChw startup__products'}).text
            else : startup_produits = ''

            #ajout dans la liste
            l_startup_nom.append(startup_nom)
            l_startup_creation.append(startup_creation)
            l_startup_identite.append(startup_identite)
            l_startup_clients.append(startup_clients)
            l_startup_produits.append(startup_produits)
            l_startup_lien.append(line)

            print("Societe : {} - C : {}".format(startup_nom, startup_creation))

        except:
            print(Fore.RED + "Erreur sur cette page." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['startup_nom'] = l_startup_nom
df_posts['startup_creation'] = l_startup_creation
df_posts['startup_identite'] = l_startup_identite
df_posts['startup_clients'] = l_startup_clients
df_posts['startup_produits'] = l_startup_produits
df_posts['startup_lien'] = l_startup_lien

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_bpi):
    os.remove(fichier_sortie_bpi)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_bpi)

#fermeture
driver.quit()