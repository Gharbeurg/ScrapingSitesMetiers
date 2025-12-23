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
fichier_entree_cardio = "C:/DATA/github/.params/liens_cardio.txt"
fichier_sortie_cardio = "C:/DATA/github/.data/cardio_content_detail.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
temporisation = 4

df_posts = pd.DataFrame(columns =  ['article_nom', 'article_tag', 'article_lien'])
df_posts = df_posts.reset_index(drop=True)
l_article_lien = []
l_article_nom = []
l_article_tag = []

# OUverture d'un VPN
print("{} - Ouverture du VPN".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
settings = initialize_VPN(area_input=['Belgium,France,Netherlands,Germany,Slovenia,Switzerland,Denmark,Croatia,Finland,Estonia,Portugal,Greece, Serbia,Austria,Spain,Poland'])
rotate_VPN(settings) #refer to the instructions variable here
time.sleep(temporisation)

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# parsing de chaque startup
print("{} - Parsing de chaque article".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

with open(fichier_entree_cardio, 'r', encoding="utf-8") as f:
    for line in f:
        nombre_requete +=1

        driver.get(line)
        time.sleep(temporisation)

        #Recuperation du contenu
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

        try:
            # Nom de la article
            if soup_contenu.find('h1', attrs = {'class': 'page-title'}):
                article_nom = soup_contenu.find('h1', attrs = {'class': 'page-title'}).text
            else : article_nom = ''
    
            # tag de l'article
            if soup_contenu.find('div', attrs = {'class': 'page-date-and-tags'}):
                article_tag = soup_contenu.find('div', attrs = {'class': 'page-date-and-tags'}).text
            else : article_tag = ''

            #ajout dans la liste
            l_article_nom.append(article_nom)
            l_article_tag.append(article_tag)
            l_article_lien.append(line)

            print("Article : {}".format(article_nom))

        except:
            print(Fore.RED + "Erreur sur cette page." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['article_nom'] = l_article_nom
df_posts['article_tag'] = l_article_tag
df_posts['article_lien'] = l_article_lien

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_cardio):
    os.remove(fichier_sortie_cardio)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_cardio)

#fermeture
driver.quit()