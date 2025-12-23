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
url_mcarthur = "https://www.ellenmacarthurfoundation.org/resources/business/circular-startup-index?csi-index[page]="
fichier_sortie_mcarthur = "C:/DATA/github/.data/mcarthur_content.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 81
temporisation = 4

df_posts = pd.DataFrame(columns =  ['name', 'country', 'text', 'tags'])
df_posts = df_posts.reset_index(drop=True)
l_post_name = []
l_post_country = []
l_post_text = []
l_post_tags = []

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# ouverture du driver
driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# Collecte des sociétés
print("{} - Démarrage de la collecte".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

for i in range(1, page_max):
    
    # parsing de chaque compte
    print("{} - page {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i))

    # Construction de la requête de recherche
    adresse = url_mcarthur + str(i)
    driver.get(adresse)
    time.sleep(temporisation)

    #Recuperation du contenu
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    # societes
    liste_societes = soup_contenu.find_all("div", {"class":"css-158rhkc e17ci1cf1"})
    for societe in liste_societes:
        try:
            # Nom
            if societe.find('p', attrs = {'class': 'startup-name'}):
                societe_nom = societe.find('p', attrs = {'class': 'startup-name'}).text
            else : societe_nom = ''

            # pays
            if societe.find('span', attrs = {'class': 'tile-country'}):
                societe_pays = societe.find('span', attrs = {'class': 'tile-country'}).text
            else : societe_pays = ''

            # texte
            if societe.find('p', attrs = {'class': 'startup-text'}):
                societe_texte = societe.find('p', attrs = {'class': 'startup-text'}).text
            else : societe_texte = ''
            
            # tags
            if societe.find('div', attrs = {'class': 'tile-work-areas'}):
                societe_tags = societe.find('div', attrs = {'class': 'tile-work-areas'}).text
            else : societe_tags = ''

            #ajout dans la liste
            l_post_name.append(societe_nom)
            l_post_country.append(societe_pays)
            l_post_text.append(societe_texte)
            l_post_tags.append(societe_tags)

        except:
            print(Fore.RED + "Erreur sur ce post." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['name'] = l_post_name
df_posts['country'] = l_post_country
df_posts['text'] = l_post_text
df_posts['tags'] = l_post_tags

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_mcarthur):
    os.remove(fichier_sortie_mcarthur)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_mcarthur)

#fermeture
driver.quit()