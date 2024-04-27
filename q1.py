import numpy as np
import matplotlib.pyplot as plt

# Constantes
R = 8.314  # Constante universelle des gaz parfaits en J/(mol·K)
rho_s = 2550  # Densité du lit en kg/m^3
T = 340 + 273.15  # Convert temperature from °C to K
P_atm = 25  # Pressure in atm
P = P_atm * 101325  # Convert pressure to Pa
dp = 0.005
d = 3  # Reactor diameter in m
L = 10  # Reactor length in m
n_sections = 100  # Nombre de sections dans le réacteur
dz = L / n_sections  # Longueur de chaque section
M_CO = 12 + 16
M_CO2 = 12 + 16 + 16
M_N2 = 14 +14
M_H2 = 2
M_H2O = 2 + 16
M_tot = M_CO+M_CO2+M_H2+M_H2O+M_N2
e_void = 0.38 + 0.073*(1+(((d/dp)-2)**2/(d/dp)**2))
A_reactor = np.pi * (d/2)**2
rho_bed = (1-e_void)*rho_s
cp = 2800
deltaH = - 38000

# Flux d'entrées pour chaque espèce en mol/s
F_CO = 240
F_CO2 = 60
F_H2 = 780
F_H2O = 410
F_N2 = 340


#=======================================================
#========================= Q1 ==========================
#=======================================================


# Total molar flux
F_tot = F_CO + F_CO2 + F_N2 + F_H2 + F_H2O

# Calculate volumetric flow rate (Q) using ideal gas law [m³/s]
Q = (F_tot * R * T) / P
u_s = Q/A_reactor
print("u_s : ",u_s)

# Calculate initial concentrations in mol/m^3
C_CO0 = F_CO / Q
C_CO20 = F_CO2 / Q
C_N20 = F_N2 / Q
C_H20 = F_H2 / Q
C_H2O0 = F_H2O / Q

# Concentrations initiales
concentrations = {
    'CO': [C_CO0] + [0] * n_sections,
    'CO2': [C_CO20] + [0] * n_sections,
    'H2': [C_H20] + [0] * n_sections,
    'H2O': [C_H2O0] + [0] * n_sections,
    'N2': [C_N20] + [0] * n_sections
}

conversion = [0] * (n_sections+1)

# Définir le taux de réaction (rate_reaction)
def rate_reaction(C_CO, C_H2O, C_CO2, C_H2, T):
    k = 7.5 * np.exp((-96000/R) * ((1/T) - (1/673.15)))  # Coefficient cinétique
    K = 10 ** ((1910/T) - 1.764)  # Constante d'équilibre
    K_CO = 0.6
    K_H2O = 1.3
    K_CO2 = 4.1
    K_H2 = 4.1
    # Calcul du taux de réaction
    numerator = k * K_CO * K_H2O * (C_CO * C_H2O - C_CO2 * C_H2 / K)
    denominator = (1 + K_CO * C_CO + K_H2O * C_H2O + K_CO2 * C_CO2 + K_H2 * C_H2) ** 2
    return numerator / denominator

# Boucle pour parcourir chaque section du réacteur

for i in range(n_sections):
    # Récupérer les concentrations actuelles
    C_CO = concentrations['CO'][i]
    C_CO2 = concentrations['CO2'][i]
    C_H2 = concentrations['H2'][i]
    C_H2O = concentrations['H2O'][i]
    C_N2 = concentrations['N2'][i]

    # Calculer le taux de réaction
    r = rate_reaction(C_CO, C_H2O, C_CO2, C_H2, T)

    # Calculer les variations de concentration pour chaque espèce
    dC_CO = 0.5* (-dz) * rho_bed * r / u_s
    dC_CO2 = 0.5* dz * rho_bed * r / u_s
    dC_H2 = 0.5 * dz * rho_bed * r / u_s
    dC_H2O = 0.5 * (-dz) * rho_bed * r / u_s
    dC_N2 = 0  # N2 n'est pas réactif, donc pas de variation de concentration

    # Mettre à jour les concentrations des espèces dans la section suivante
    concentrations['CO'][i + 1] = C_CO + dC_CO
    concentrations['CO2'][i + 1] = C_CO2 + dC_CO2
    concentrations['H2'][i + 1] = C_H2 + dC_H2
    concentrations['H2O'][i + 1] = C_H2O + dC_H2O
    concentrations['N2'][i + 1] = C_N2
    conversion[i] = ((concentrations['CO'][0]-concentrations['CO'][i])/concentrations['CO'][0])


# Afficher les concentrations des espèces en fonction de la longueur du réacteur
longueur = np.linspace(0, 10, n_sections + 1)
plt.plot(longueur, concentrations['CO'], label='CO')
plt.plot(longueur, concentrations['CO2'], label='CO2')
plt.plot(longueur, concentrations['H2'], label='H2')
plt.plot(longueur, concentrations['H2O'], label='H2O')
plt.plot(longueur, concentrations['N2'], label='N2')

plt.xlabel('Lenght of reactor (m)')
plt.ylabel('Concentration ($mol/m^3$)')
plt.title('Concentrations of species along the reactor')
plt.legend()
plt.show()
conversion[-1] = ((concentrations['CO'][0]-concentrations['CO'][-1])/concentrations['CO'][0])
conversion_CO = ((concentrations['CO'][0]-concentrations['CO'][-1])/concentrations['CO'][0])
print("Conversion rate of CO = ",conversion_CO*100,"%")
print ("Void fraction = ",e_void*100,"%")

plt.plot(longueur,conversion,label = 'CO conversion')
plt.xlabel('Lenght of reactor (m)')
plt.ylabel('CO conversion (-)')
plt.title('Conversion of CO along the reactor')
plt.legend()
plt.show()
#=======================================================
#========== Second way to plot results of Q1 ===========
#=======================================================
n_tot = concentrations['CO'][0]+concentrations['CO2'][0]+concentrations['N2'][0]+concentrations['H2'][0]+concentrations['H2O'][0]
for i in range(len(concentrations['CO'])):
    concentrations['CO'][i] = concentrations['CO'][i]/n_tot
    concentrations['CO2'][i] = concentrations['CO2'][i]/n_tot
    concentrations['N2'][i] = concentrations['N2'][i]/n_tot
    concentrations['H2'][i] = concentrations['H2'][i]/n_tot
    concentrations['H2O'][i] = concentrations['H2O'][i]/n_tot
plt.plot(longueur, concentrations['CO'], label='CO')
plt.plot(longueur, concentrations['CO2'], label='CO2')
plt.plot(longueur, concentrations['H2'], label='H2')
plt.plot(longueur, concentrations['H2O'], label='H2O')
plt.plot(longueur, concentrations['N2'], label='N2')

plt.xlabel('Length of reactor')
plt.ylabel('Molar fraction (-)')
plt.title('Molar fraction along the reactor')
plt.legend()
plt.show()