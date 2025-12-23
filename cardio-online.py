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
url_cardio = "https://www.cardio-online.fr/Actualites/(offset)/"
fichier_sortie_cardio = "C:/DATA/github/.data/cardio_content.xlsx"
nbre_pages_erreur = 0
nombre_requete = 0
page_max = 30
temporisation = 4

df_articles = pd.DataFrame(columns =  ['article_titre', 'article_date', 'article_theme', 'article_auteur', 'article_lien'])
df_articles = df_articles.reset_index(drop=True)
l_article_titre = []
l_article_date = []
l_article_theme = []
l_article_auteur = []
l_article_lien = []

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
print("{} - Démarrage de la collecte des articles".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

for i in range(0, page_max):
    
    # parsing de chaque compte
    print("{} - page {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i))

    # Construction de la requête de recherche
    adresse = url_cardio + str(i) + "0"
    driver.get(adresse)
    time.sleep(temporisation)

    #Recuperation du contenu
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    # liste des liens
    liste_articles = soup_contenu.find_all("div", {"class":"media embed-listitem"})
    for article in liste_articles:

        try:
            # données de l'article
            if article.find("h4", {"class":"media-heading"}):
                titre_article = article.find("h4", {"class":"media-heading"}).text
            else: titre_article = ''

            if article.find("a"):
                lien_article = article.find("a").get('href')
                lien_article = "https://www.cardio-online.fr" + lien_article
            else: lien_article = ''
            
            if article.find("p", {"class":"publication-date"}):
                date_article = article.find("p", {"class":"publication-date"}).text
            else: date_article = ''
            
            theme_article = ''
            
            if article.find("div", {"class":"media-description"}):
                auteur_article = article.find("div", {"class":"media-description"}).text
            else: auteur_article = ''
            
            #ajout dans la liste
            l_article_titre.append(titre_article)
            l_article_date.append(date_article)
            l_article_theme.append(theme_article)
            l_article_auteur.append(auteur_article)
            l_article_lien.append(lien_article)

        except:
            print(Fore.RED + "Erreur sur cette page." + Fore.RESET)
            
# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles['article_titre'] = l_article_titre
df_articles['article_date'] = l_article_date
df_articles['article_theme'] = l_article_theme
df_articles['article_auteur'] = l_article_auteur
df_articles['article_lien'] = l_article_lien

# ecriture des fichiers de sortie
if os.path.exists(fichier_sortie_cardio):
    os.remove(fichier_sortie_cardio)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles.to_excel(fichier_sortie_cardio)

#fermeture
driver.quit()