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
url_quotidien_medecin = "https://www.lequotidiendumedecin.fr/actus-medicales/esante?page="
fichier_sortie_quotidien = "C:/DATA/github/.data/quotidien_content.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 6
temporisation = 4

df_posts = pd.DataFrame(columns =  ['post_titre', 'post_date', 'post_description', 'post_nombre_reaction', 'post_lien'])
df_posts = df_posts.reset_index(drop=True)
l_post_titre = []
l_post_date = []
l_post_description = []
l_post_reaction = []
l_post_lien = []

# OUverture d'un VPN
print("{} - Ouverture du VPN".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
settings = initialize_VPN(area_input=['Belgium,France,Netherlands,Germany,Slovenia,Switzerland,Denmark,Croatia,Finland,Estonia,Portugal,Greece, Serbia,Austria,Spain,Poland'])
rotate_VPN(settings) #refer to the instructions variable here
time.sleep(temporisation)

# Ouverture du browser
print("{} - Ouverture du driver".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# ouverture du driver
driver = webdriver.Chrome(executable_path=r'C:/DATA/github/drivers/chromedriver.exe')

# Collecte des posts
print("{} - Démarrage de la collecte".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

for i in range(0, page_max):
    
    # parsing de chaque compte
    print("{} - page {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i))

    # Construction de la requête de recherche
    adresse = url_quotidien_medecin + str(i)
    driver.get(adresse)
    time.sleep(temporisation)

    #Recuperation du contenu
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    # articles
    liste_posts = soup_contenu.find_all("div", {"class":"layout layout--gps-listing"})
    for post in liste_posts:
        try:
            # titre post
            if post.find('a'):
                titre_post = lien_post = post.find('a').get('title')
            else : titre_post = ''

            # date du post
            if post.find('div', attrs = {'class': 'font-size-16 text-uppercase text-gray-500 d-none d-md-inline'}):
                date_post = post.find('div', attrs = {'class': 'font-size-16 text-uppercase text-gray-500 d-none d-md-inline'}).text
                date_post = date_post[:-3]
            else : date_post = ''


            # description du post
            if post.find('span', attrs = {'class': 'text'}):
                description_post = post.find('span', attrs = {'class': 'text'}).text
                description_post = description_post.strip() + '\n'
            else : description_post = ''

            # réactions du post
            if post.find('span', attrs = {'class': 'nb-com'}):
                reaction_post = post.find('span', attrs = {'class': 'nb-com'}).text
                
                #suppression des espaces
                reaction_post = re.sub(' +', ' ', reaction_post)
                reaction_post = reaction_post.strip() + '\n'
            else : reaction_post = ''

            # lien du post
            if post.find('a'):
                lien_post = post.find('a').get('href')
                lien_post = 'https://www.lequotidiendumedecin.fr' + lien_post
            else : lien_post = ''

            #ajout dans la liste
            l_post_titre.append(titre_post)
            l_post_date.append(date_post)
            l_post_description.append(description_post)
            l_post_reaction.append(reaction_post)
            l_post_lien.append(lien_post)

        except:
            print(Fore.RED + "Erreur sur ce post." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts['post_titre'] = l_post_titre
df_posts['post_date'] = l_post_date
df_posts['post_description'] = l_post_description
df_posts['post_nombre_reaction'] = l_post_reaction
df_posts['post_lien'] = l_post_lien

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_quotidien):
    os.remove(fichier_sortie_quotidien)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_posts.to_excel(fichier_sortie_quotidien)

#fermeture
driver.quit()