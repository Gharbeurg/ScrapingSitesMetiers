# Bibliotheque
import os
import pandas as pd
import re
import bs4

from datetime import datetime

pattern = re.compile(r'\s+')
fichier_entree = "C:/DATA/github/.data/pubmed_gaucher.html"
fichier_pubmed_sortie = "C:/DATA/github/.data/pubmed_sortie.xlsx"
df_pages = pd.DataFrame(columns =  ['id', 'if', 'quartile', 'titre', 'auteurs', 'journal', 'creation', 'keywords', 'resume'])
Mois_publication = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
l_id = []
l_if = []
l_quartile = []
l_titre = []
l_auteurs = []
l_journal = []
l_creation = []
l_background =[]
l_aim = []
l_methods = []
l_results = []
l_conclusion = []
l_keywords = [] 
l_resume =[]
token = 0

#Ouverture du fichier d'entrée
print("{} - Ouverture du fichier d'entrée".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f = open(fichier_entree, 'r', encoding="utf-8")
soup_contenu = bs4.BeautifulSoup(f, 'html.parser')

#Liste des articles
print("{} - Liste des articles".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
tag = soup_contenu.find_all("div", {"class":"docsum-wrap"})

#Parsing de chaque article
for el in tag:
  
    #ID article
    tag = el.find("span", {"class":"docsum-pmid"})
    if tag is not None:
        print("{} - Parsing article : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), tag.text))
        l_id.append(tag.text)
    else:
        l_id.append("")

    #IF article
    tag = el.find("span", {"class":"if"})
    if tag is not None:
        l_if.append(tag.text.replace(".", ","))
    else:
        l_if.append("")

    #quartile article
    tag = el.find("span", {"class":"quartile"})
    if tag is not None:
        l_quartile.append(tag.text)
    else:
        l_quartile.append("")    

    #Titre article
    tag = el.find("a", {"class":"docsum-title"})
    if tag is not None:
        l_titre.append(tag.text.strip())
    else:
        l_titre.append("")    

    #auteurs
    tag = el.find("span", {"class":"docsum-authors full-authors"})
    if tag is not None:
        #suppression des caractères en trop
        liste_auteurs = ''.join(i for i in tag.text if not i.isdigit())
        liste_auteurs = " ".join(liste_auteurs.split())
        liste_auteurs = liste_auteurs.replace("\n", "")
        l_auteurs.append(liste_auteurs)
    else:
        l_auteurs.append("")    

    #journal
    tag = el.find("i")
    if tag is not None:
         l_journal.append(tag.text.strip())
    else:
        l_journal.append("")    

    #creation
    tag = el.find("span", {"class":"docsum-journal-citation full-journal-citation"})
    if tag is not None:
        index_date = tag.text.find(').')
        date_creation = tag.text[index_date + 3 :index_date + 11]
        l_creation.append(date_creation)
    else:
        l_creation.append("")    

    #resume
    tag = el.find("div", {"class":"full-view-snippet"})
    if tag is not None:
        # resume
        text_resume = tag.text.strip()
        text_resume = " ".join(text_resume.split())
        text_resume = text_resume.replace("\n", "")
        l_resume.append(text_resume)
    else:
        l_resume.append("")   

    #keywords
    tag = el.find("div", {"class":"short-view-snippet"})
    if tag is not None:
        # resume
        text_keywords = tag.text.strip()
        text_keywords = " ".join(text_keywords.split())
        text_keywords = text_keywords.replace("\n", "")
        l_keywords.append(text_keywords)
    else:
        l_keywords.append("")   

#creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_pages['id'] = l_id
df_pages['if'] = pd.Series(l_if)
df_pages['quartile'] = pd.Series(l_quartile)
df_pages['titre'] = pd.Series(l_titre)
df_pages['auteurs'] = pd.Series(l_auteurs)
df_pages['journal'] = pd.Series(l_journal)
df_pages['creation'] = pd.Series(l_creation)
df_pages['keywords'] = pd.Series(l_keywords) 
df_pages['resume'] = pd.Series(l_resume)

# ecriture des fichiers de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_pubmed_sortie):
    os.remove(fichier_pubmed_sortie)

df_pages.to_excel(fichier_pubmed_sortie)