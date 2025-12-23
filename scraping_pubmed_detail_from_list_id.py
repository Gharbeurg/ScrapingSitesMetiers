#bibliotheques
import bs4
import time
import os
import unicodedata
from unidecode import unidecode
import re
from tabulate import tabulate
from datetime import datetime
from colorama import Fore, Back, Style
import pandas as pd
import cloudscraper

from nltk.tokenize import sent_tokenize

from selenium import webdriver
from selenium.webdriver.chrome.options import Options

#fichiers
fichier_entree = "C:/DATA/github/.params/entree.txt"
fichier_sortie = "C:/DATA/github/.data/content_pubmed.txt"

#variables
nbre_phrases_totales = 0
nbre_phrases_traitees = 0
nbre_phrases_erreur = 0
nbre_minimum_caractere_phrase = 20
nbre_pages = 0
nbre_pages_erreur = 0
liste_pages_web = []
doi_link = "https://doi.org/"
temporisation = 50
scraper = cloudscraper.create_scraper()   #bypasscloudflare

df_table_phrases = pd.DataFrame(columns =  ['phrase', 'page', 'fichier', 'balise'])
df_table_phrases = df_table_phrases.reset_index(drop=True)
l_phrase = []
l_fichier = []
l_page = []
l_balise = []

#driver chrome
option = webdriver.ChromeOptions()
option.add_experimental_option('excludeSwitches', ['enable-logging'])
option.add_argument('--ignore-certificate-errors')
driver = webdriver.Chrome(options=option)

# suppression du fichier de sortie
print("{} - Création du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_sortie):
    os.remove(fichier_sortie)

file_sortie = open(fichier_sortie, 'a', encoding='utf-8')

# Collecte des pages à travailler
print("{} - Lecture du fichier d'entrée".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
with open(fichier_entree, 'r') as f:
    for line in f:
        # suppression des caractères spéciaux
        line = doi_link + line.strip().lower()
        liste_pages_web.append(line)

f.close()

# parsing de chaque page web
print("{} - Récupération du contenu de chaque page".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
for i in liste_pages_web:
    try:
        #recuperation de la page
        driver.get(i)
        get_url = driver.current_url
        time.sleep(temporisation)
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')

        print("[+] page : {} ".format(get_url))

        #URL de l'article
        file_sortie.write(get_url + "\n")

        #AMEGROUP
        if "amegroups" in get_url:
            soup_article = soup_contenu.find("article",  {"class":"article"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #ASCOPUB
        if "ascopubs" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"core-container"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #BIOMED
        if "biomedcentral" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"main-content"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #CUREUS
        if "cureus" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"new-article-content-wrap"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #CELL
        if "cell" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"col-md-7 col-lg-9 article__sections"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #clinical-lung-cancer
        if "clinical-lung-cancer" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"row"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #DOVEPRESS
        if "dovepress" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"tab-content"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")


        #FRONTIERIN
        if "frontiersin" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"ArticleDetails__main__content"})
            
            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #INDIANJOURNALNEPHRO
        if "indianjnephrol" in get_url:
            soup_article = soup_contenu.find("article")

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #JTCVS
        if "jtcvs" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"row"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #JTO
        if "jto" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"container"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #LWW
        if "lww" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"ArticleContainer"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #MDPI
        if "mdpi" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"html-article-content"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #NATURE
        if "nature" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"c-article-body"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #NEJM
        if "nejm" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"core-container"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #PORJOUTNAL
        if "por-journal" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"JournalFullText"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #SAGE
        if "sagepub" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"core-container"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #SCIENCEDIRECT
        if "sciencedirect" in get_url:
            soup_article = soup_contenu.find("article",  {"class":"col-lg-12 col-md-16 pad-left pad-right"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")


        #SPANDIDOS
        if "spandidos" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"content_main"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #SPRINGER
        if "springer" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"c-article-body"})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

        #WILEY
        if "wiley" in get_url:
            soup_article = soup_contenu.find("div",  {"class":"article__body "})

            # Tokenisation
            sentence_tokens = sent_tokenize(soup_article.text.strip())

            # Nettoyage des caractères spéciaux
            for sent in sentence_tokens:
                sent = re.sub('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+#]|[!*\(\),]|''(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', sent)
                sent = re.sub("(@[A-Za-z0-9_]+)", "", sent)
                sent = sent.replace("\n", "")
                sent = sent.encode('ascii', 'ignore').decode('ascii')
                sent = sent.lower()
                sent = sent.strip()
                if len(sent) > nbre_minimum_caractere_phrase:
                    file_sortie.write(sent)
                    file_sortie.write("\n")

    except:
        print(Fore.RED + "[+] Erreur sur la page : {}".format(get_url) + Fore.RESET)
        nbre_pages_erreur += 1

print("{} - Fin du parsing des pages".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
print("[+] Nombe de pages  en erreur : ", nbre_pages_erreur)

#fermeture
driver.quit()
file_sortie.close()

print("{} - Fin du programme".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))