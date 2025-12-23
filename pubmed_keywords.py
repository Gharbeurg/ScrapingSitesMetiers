#Bibliotheques
import pandas as pd
import time
import bs4
import os

from selenium import webdriver
from nordvpn_switcher import initialize_VPN, rotate_VPN

from datetime import datetime
from colorama import Fore, Back, Style

#Variables
fichier_id = "C:/DATA/github/.data/pubmed_auteurs.txt"
token_pubmed = "https://pubmed.ncbi.nlm.nih.gov/" # adresse de pubmed
liste_motscles = ['burden', 'patient', 'patients', 'fruit', 'alcohol', 'child', 'children', 'adult', 'teenager', 'woman', 'pregnancy', 'clinical', 'therapy', 'quality', 'life', 'treatment', 'treatments', 'prophylactic', 'burden', 'emotion', 'anxiety', 'impact' 'journey', 'diagnosis', 'improving quality of life', 'travel']
fichier_pubmed_sortie = "C:/DATA/github/.data/pubmed__keywords.xlsx"

temporisation = 2 #temporisation entre chaque click en secondes
df_articles = pd.DataFrame(columns =  ['id','mot'])
df_articles = df_articles.reset_index(drop=True)
l_id = []
l_motcle = []

# Collecte des pages à travailler
print("{} - Ouverture du fichier des articles".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

# connexion à pubmed et initialisation d'un VPN
settings = initialize_VPN(area_input=['Belgium,France,Netherlands,Germany,Slovenia,Switzerland','Denmark','Croatia','Finland','Estonia'])
rotate_VPN(settings)

driver = webdriver.Chrome(executable_path="C:/DATA/github/drivers/chromedriver.exe")

#temporisation
time.sleep(temporisation)

with open(fichier_id, 'r', encoding="utf-8") as f:
    for line in f:
        try:
            # Construction de la requête
            requete = token_pubmed + line
            driver.get(requete)

            #temporisation
            time.sleep(temporisation)

            #Recuperation du contenu
            contenu = driver.page_source.encode('utf8')
            soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')
       
            result_block = soup_contenu.find_all(["h1", "h2", "p"])
            counts = dict()

            for result in result_block:
                occurence = 0
                resume = result.text.lower()
                resume_liste_mots = resume.split()

                for mot_cle in liste_motscles:
                    occurence = resume_liste_mots.count(mot_cle)
                    if occurence > 0:
                        if mot_cle in counts:
                            counts[mot_cle] += occurence
                        else:
                            counts[mot_cle] = occurence
                        
            l_id.append(line)
            l_motcle.append(counts)

            print("{} - Traitement article : {} - Mots trouvés : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), line.rstrip(), counts))
        
        except:
            print(Fore.RED + "{} - Traitement article : {} - Nada".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), line.rstrip()) + Fore.RESET)

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles['mot'] = l_motcle
df_articles['id'] = l_id

# ecriture des fichiers de sortie
if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_articles.to_excel(fichier_pubmed_sortie)

# Fermture du fichier des auteurs
print("{} - Fermeture du fichier des articles".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f.close()
driver.quit()