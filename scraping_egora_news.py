# Bibliotheque
import time
import pandas as pd
import time
import bs4
import re

from datetime import datetime
from colorama import Fore, Back, Style
from selenium import webdriver

# variables
temporisation = 5
nombre_requete = 0

adresse_driver_chrome = "C:/DATA/github/drivers/chromedriver.exe"
adresse_egora_news = "https://www.egora.fr/actus"
adresse_egora = "https://www.egora.fr"
fichier_egora_sortie = "C:/DATA/github/.data/sante_news.csv"
fichier_log = "C:/DATA/github/.logs/log_execution.txt"
journal = 'EGORA'

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
f_log.write("{} - Lancement du programme - scraping_egora_news.py\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# Ouverture du Browser
f_log.write("{} - Ouverture du browser\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

options = webdriver.ChromeOptions() 
options.add_experimental_option("excludeSwitches", ["enable-logging"])
driver = webdriver.Chrome(options=options, executable_path=adresse_driver_chrome)

# ouverture de la page
f_log.write("{} - Ouverture de la page\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
try:
    driver.get(adresse_egora_news)
    time.sleep(temporisation)

    #Recuperation du contenu
    f_log.write("{} - Recuperation du contenu\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    contenu = driver.page_source.encode('utf8')
    soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

    #récuperation des articles
    regex = re.compile('.*views-row views-row-.*')
    result_block = soup_contenu.find_all("div", {"class":regex})

    if result_block:
        try:
            for result in result_block:
                nombre_requete += 1
        
                #journal_article
                journal_article = journal

                #categorie
                if result.find("span", {"class":"rubrique_liste"}):
                    categorie_article = result.find("span", {"class":"rubrique_liste"}).text.strip()

                else:
                    categorie_article = ''

                #titre
                if result.find("h2"):
                    titre_article = result.find("h2").text.strip()
                    commente_article = titre_article[-3:].strip()
                    titre_article = titre_article[:-3]
                else:
                    titre_article = ''
                    commente_article = ''

                #description
                if result.find("div", {"class":"field field-name-body field-type-text-with-summary field-label-hidden"}):
                    description_article = result.find("div", {"class":"field field-name-body field-type-text-with-summary field-label-hidden"}).text.strip()
                    description_article = description_article[:-10]
                else:
                    description_article = ''

                #lien
                if result.find("a", {"class":"more-link"}):
                    lien_article = adresse_egora + result.find("a", {"class":"more-link"}).get('href')
                else:
                    lien_article = ''

                #date
                if result.find("span", {"class":"date_liste"}):
                    date_article = result.find("span", {"date_liste"}).text.strip()
                else:
                    date_article = ''

                #partage
                partage_article = ''

                #auteur
                auteur_article = ''

                #ajout dans la liste
                l_journal.append(journal_article)
                l_date.append(date_article)
                l_categorie.append(categorie_article)
                l_titre.append(titre_article)
                l_auteur.append(auteur_article)
                l_description.append(description_article)
                l_lien.append(lien_article)
                l_commente.append(commente_article)
                l_partage.append(partage_article)

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
df_news.to_csv(fichier_egora_sortie, mode='a', sep=';', index=False, header=False)

#fermeture
driver.quit()
f_log.close

#fin du programme
f_log.write("{} - Fin du programme - scraping_egora_news.py\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f_log.write("****************************************************\n")