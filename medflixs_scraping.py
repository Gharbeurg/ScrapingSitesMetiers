# Bibliotheque
import os
import pandas as pd
import re
import bs4

from datetime import datetime

# Variables
pattern = re.compile(r'\s+')
fichier_entree = "C:/DATA/github/.params/agenda_congres.htm"
fichier_event_sortie = "C:/DATA/github/.data/agenda_sortie.xlsx"
df_pages = pd.DataFrame(columns =  ['startdate', 'enddate', 'nom', 'specialite', 'localisation', 'lien_event', 'lien_detail'])
l_startdate = []
l_enddate = []
l_nom = []
l_specialite = []
l_localisation = []
l_lien_event = []
l_lien_detail = []

nom_evenement = ""
i = 0

#Ouverture du fichier d'entrée
print("{} - Ouverture du fichier d'entrée".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
f = open(fichier_entree, 'r', encoding="utf-8")
soup_contenu = bs4.BeautifulSoup(f, 'html.parser')

#Liste des évènements
print("{} - collecte de la Liste des évènements".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
tag = soup_contenu.find_all("div", {"class":"event__item"})

#Parsing de chaque evenement
for el in tag:
  
    #compteur
    i += 1

    #start date
    tag = el.find("time", {"itemprop":"startDate"})
    if tag is not None:
        l_startdate.append(tag.text)
    else:
        l_startdate.append("")

    #end date
    tag = el.find("time", {"itemprop":"endDate"})
    if tag is not None:
        l_enddate.append(tag.text)
    else:
        l_enddate.append("")

    #Lien evenement
    tag = el.find("a", {"class":"link"})
    if tag is not None:
        l_lien_event.append(tag.get('href'))
    else:
        l_lien_event.append("")    

    #Lien detail
    tag = el.find("a", {"class":"event__href"})
    if tag is not None:
        l_lien_detail.append(tag.get('href'))
    else:
        l_lien_detail.append("")    

    #Nom evenement
    tag = el.find("p", {"itemprop":"name"})
    if tag is not None:
        l_nom.append(tag.text.strip())
        nom_evenement = tag.text.strip()
    else:
        l_nom.append("")
        nom_evenement = ""

    print("[{}] evenement : {}".format(i, nom_evenement))

    #specialite evenement
    tag = el.find("p", {"class":"event__spe"})
    if tag is not None:
        specialite_sans_saut_ligne = tag.text.replace("\n", "")
        specialite_sans_saut_ligne = specialite_sans_saut_ligne.replace("  ", "")
        specialite_sans_saut_ligne = specialite_sans_saut_ligne[:len(specialite_sans_saut_ligne)-1]
        l_specialite.append(specialite_sans_saut_ligne)
    else:
        l_specialite.append("")

    #localisation evenement
    tag = el.find("span", {"itemprop":"address"})
    if tag is not None:
        localisation_sans_saut_ligne = tag.text.replace("\n", "")
        localisation_sans_saut_ligne = localisation_sans_saut_ligne.replace("  ", "")
        l_localisation.append(localisation_sans_saut_ligne)
    else:
        l_localisation.append("")

#creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_pages['startdate'] = l_startdate
df_pages['enddate'] = pd.Series(l_enddate)
df_pages['nom'] = pd.Series(l_nom)
df_pages['specialite'] = pd.Series(l_specialite)
df_pages['localisation'] = pd.Series(l_localisation)
df_pages['lien_event'] = pd.Series(l_lien_event)
df_pages['lien_detail'] = pd.Series(l_lien_detail)

# ecriture des fichiers de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
if os.path.exists(fichier_event_sortie):
    os.remove(fichier_event_sortie)

df_pages.to_excel(fichier_event_sortie)