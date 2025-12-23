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
nbre_pages_erreur = 0
nbre_labos_erreur = 0
compteur = 0

adresse_driver_chrome = "C:/DATA/github/drivers/chromedriver.exe"
leem_annuaire = "https://www.leem.org/annuaire-des-entreprises-adherentes?activite=All&ordre=All&recherche=&page="
fichier_leem_sortie = "C:/DATA/github/.data/leem_annuaire.xlsx"
nbre_pages_annuaire_leem = 46

df_leem = pd.DataFrame(columns =  ['nom', 'societe', 'site', 'rue', 'postal', 'ville', 'rueplus', 'phone', 'mail', 'web'])
df_leem = df_leem.reset_index(drop=True)
l_nom = []
l_societe = []
l_site = []
l_rue = []
l_postal = []
l_ville = []
l_rueplus = []
l_phone = []
l_mail = []
l_web = []

# Construction de la liste de pages
def get_pages(token, nb):
    pages = []
    for i in range(0,nb):
        j = token + str(i)
        pages.append(j)
    return pages

print("{} - Construction de la liste de pages".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
pages = get_pages(leem_annuaire,nbre_pages_annuaire_leem)

options = webdriver.ChromeOptions() 
options.add_experimental_option("excludeSwitches", ["enable-logging"])
driver = webdriver.Chrome(options=options, executable_path=adresse_driver_chrome)

# parsing de chaque page sujet
print("{} - Collecte des informations de chaque page".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
for i in pages:
    try:
        compteur += 1
        print("{} - Traitement de la page : {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), compteur))
        driver.get(i)
        contenu = driver.page_source.encode('utf8')
        soup_contenu = bs4.BeautifulSoup(contenu, 'html.parser')
        time.sleep(temporisation)

        labos = soup_contenu.find_all("div", {"class":"views-row"})
        for labo in labos:
            if labo.find("div", {"class":"views-field views-field-title"}):
                titre = labo.find("div", {"class":"views-field views-field-title"}).text
                print(titre)
            else:
                 titre = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-entreprise"}):
                societe = labo.find("div", {"class":"views-field views-field-field-leem-entreprise"}).text
            else:
                societe = ''            
            
            if labo.find("div", {"class":"views-field views-field-field-leem-activity-type"}):
                site = labo.find("div", {"class":"views-field views-field-field-leem-activity-type"}).text
            else:
                site = ''            
            
            if labo.find("div", {"class":"views-field views-field-field-leem-street"}):
                rue = labo.find("div", {"class":"views-field views-field-field-leem-street"}).text
            else:
                rue = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-postal-code"}):    
                postal = labo.find("div", {"class":"views-field views-field-field-leem-postal-code"}).text
            else:
                postal = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-city"}):    
                ville = labo.find("div", {"class":"views-field views-field-field-leem-city"}).text
            else:
                ville = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-address-plus"}):    
                rueplus = labo.find("div", {"class":"views-field views-field-field-leem-address-plus"}).text
            else:
                rueplus = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-phone"}):
                phone = labo.find("div", {"class":"views-field views-field-field-leem-phone"}).text
            else:
                phone = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-email"}):    
                mail = labo.find("div", {"class":"views-field views-field-field-leem-email"}).text
            else:
                mail = ''
            
            if labo.find("div", {"class":"views-field views-field-field-leem-website"}):    
                web = labo.find("div", {"class":"views-field views-field-field-leem-website"}).text
            else:
                web = ''
            
            l_nom.append(titre)
            l_societe.append(societe)
            l_site.append(site)
            l_rue.append(rue)
            l_postal.append(postal)
            l_ville.append(ville)
            l_rueplus.append(rueplus)
            l_phone.append(phone)
            l_mail.append(mail)
            l_web.append(web)

            print("[+] labo : ", titre)

    except:
        nbre_labos_erreur += 1

print("[+] Nombe de labos en erreur sur cette page : ", nbre_labos_erreur )
nbre_labos_erreur = 0

# creation du dataframe
print("{} - Création du dataframe".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
df_leem['nom'] = l_nom
df_leem['societe'] = l_societe
df_leem['site'] = l_site
df_leem['rue'] = l_rue
df_leem['postal'] = l_postal
df_leem['ville'] = l_ville
df_leem['rueplus'] = l_rueplus
df_leem['phone'] = l_phone
df_leem['mail'] = l_mail
df_leem['web'] = l_web


# ecriture des fichiers de sortie
print("{} - Ecriture du fichier de sortie".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
with open(fichier_leem_sortie, 'w', encoding="utf-8") as fs:
    df_leem.to_excel(fichier_leem_sortie)

#fermeture
driver.quit()

#fin du programme
print("{} - Fin du programme".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))