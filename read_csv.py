import pandas as pd
import matplotlib.pyplot as plt

# Chemin vers ton fichier Excel
chemin_fichier_excel = "/Users/hugosimon/Desktop/PROJET REACTOR DESIGN/Fixed_bed/fanning coeff.csv.xlsx"

# Lire le fichier Excel
donnees_excel = pd.read_excel(chemin_fichier_excel)

# Récupérer les données des trois colonnes dans des listes Python
length = donnees_excel.iloc[:, 0].tolist()
fanning1 = donnees_excel.iloc[:, 1].tolist()
fanning2 = donnees_excel.iloc[:, 2].tolist()
fanning3 = donnees_excel.iloc[:, 3].tolist()

# Afficher les trois premiers éléments de chaque colonne pour vérification
print("Première colonne :", length[:3])
print("Deuxième colonne :", fanning1[:3])
print("Troisième colonne :", fanning2[:3])
print("Quatrième colonne :", fanning3[:3])

plt.plot(length, fanning1, label='Ergun + Burke and Plummer')
plt.plot(length, fanning2, label='Ergun')
plt.plot(length, fanning3, label='Handley and Heggs')


# Configurer le graphique
plt.xlabel('Lenght of reactor (m)')
plt.ylabel('Pressure [atm]')
plt.title('Pressure along the reactor')
plt.legend()
plt.show()