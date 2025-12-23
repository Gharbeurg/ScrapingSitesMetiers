# Bibliotheque
import time
import pandas as pd
import time
import bs4

from datetime import datetime
from colorama import Fore, Back, Style
from selenium import webdriver

# variables
temporisation = 5
nombre_requete = 0

adresse_driver_chrome = "C:/DATA/github/drivers/chromedriver.exe"
adresse_univadis_news = "https://www.univadis.fr/news/all/"
adresse_univadis = "https://www.univadis.fr"
fichier_univadis_sortie = "C:/DATA/github/.data/sante_news.csv"
fichier_log = "C:/DATA/github/.logs/log_execution.txt"
journal = 'UNIVADIS-FR'

df_news = pd.DataFrame(columns =  ['journal', 'date', 'categorie', 'auteur', 'titre', 'description', 'lien', 'commente', 'partage'])
df_news = df_news.reset_index(drop=True)
l_journal = []
l_date = []
l_categorie = []
l_titre = []
l_auteur = []
l_description = []
l_lien = []
l_commente = []
l_partage = []

# ouverture du fichier de logs
f_log = open(fichier_log, 'a', encoding='utf8')
f_log.write("{} - Lancement du programme - scraping_univadis_news.py\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# Ouverture du Browser
f_log.write("{} - Ouverture du browser\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

options = webdriver.ChromeOptions() 
options.add_experimental_option("excludeSwitches", ["enable-logging"])
driver = webdriver.Chrome(options=options, executable_path=adresse_driver_chrome)

# ouverture de la page
f_log.write("{} - Ouverture de la page\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
try:
    driver.get(adresse_univadis_news)
    time.sleep(temporisation)

    #Recuperation du contenu
    f_log.write("{} - Recuperation du contenu\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    #récuperation des articles
    result_block = soup_contenu.find_all("li", {"class":"new-homepage-article__item"})

    if result_block:
        try:
            for result in result_block:
                nombre_requete += 1
        
                #journal_article
                journal_article = journal

                #meta, auteur, categorie
                if result.find("span", {"class":"new-homepage-article__meta"}):
                    meta_article = result.find("span", {"class":"new-homepage-article__meta"}).text.rstrip()
                    limite = meta_article.find('-')
                    date_article = meta_article[:limite-1].strip()
                    categorie_article = meta_article[limite+1:].strip()
                    limite = categorie_article.find(' de ')
                    auteur_article = categorie_article[limite+3:].strip()
                    categorie_article = categorie_article[:limite].strip()

                else:
                    meta_article = ''
                    date_article = ''
                    categorie_article = ''
                    auteur_article = ''

                #titre
                if result.find("strong", {"class":"new-homepage-article__subtitle"}):
                    titre_article = result.find("strong", {"class":"new-homepage-article__subtitle"}).text.rstrip()
                else:
                    titre_article = ''

                #description
                if result.find("p", {"class":"new-homepage-article__short-line"}):
                    description_article = result.find("p", {"class":"new-homepage-article__short-line"}).text.rstrip()
                else:
                    description_article = ''

                #lien
                if result.find("a", {"class":"new-homepage-article__article-link"}):
                    lien_article = adresse_univadis + result.find("a", {"class":"new-homepage-article__article-link"}).get('href')
                else:
                    lien_article = ''

                #ajout dans la liste
                l_journal.append(journal_article)
                l_date.append(date_article)
                l_categorie.append(categorie_article)
                l_titre.append(titre_article)
                l_auteur.append(auteur_article)
                l_description.append(description_article)
                l_lien.append(lien_article)
                l_commente.append("")
                l_partage.append("")

                f_log.write("[+] - {} - {}\n".format(nombre_requete, titre_article.rstrip()))
        
        except:
            f_log.write(Fore.RED + "[+] - Erreur docteur, ligne {}\n".format(nombre_requete) + Fore.RESET)

except:
    f_log.write(Fore.RED + "[+] - Mince ca deconne\n" + Fore.RESET)

# creation du dataframe
f_log.write("{} - Création du dataframe\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_news['journal'] = l_journal
df_news['date'] = l_date
df_news['categorie'] = l_categorie
df_news['titre'] = l_titre
df_news['auteur'] = l_auteur
df_news['description'] = l_description
df_news['lien'] = l_lien
df_news['commente'] = l_commente
df_news['partage'] = l_partage

# ecriture des fichiers de sortie
f_log.write("{} - Ecriture du fichier de sortie\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_news.to_csv(fichier_univadis_sortie, mode='a', sep=';', index=False, header=False)

#fermeture
driver.quit()
f_log.close

#fin du programme
f_log.write("{} - Fin du programme - scraping_univadis_news.py\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f_log.write("****************************************************\n")