import math
import csv
from scipy.optimize import fsolve
import numpy as np
from scipy.constants import G
from scipy.constants import Stefan_Boltzmann as sigma


#berechnung von x_o

def height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o):

    a = -(p_o-p_u)
    b = 0
    c = -3*g_p*M_H2O/(4*math.pi)
    d = 3*g_p*M_H2O*x_u/(4*math.pi) -a*x_u**3

    # Koeffizienten in einer Liste, um sie an np.roots zu übergeben
    coefficients = [a, b, c, d]

    # Berechnen der Wurzeln
    roots = np.roots(coefficients)
    print(roots)
    
    hl = []
    
    for i in range(0,3):
        if roots[i]>x_u:
            hl.append(roots[i])
    #print(hl)
    x_o = min(hl)
    
    return x_o

#berechnung von x_m

def height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o):
    x_m = (3*M_H2O/(4*math.pi*(rho_w-rho_e))-rho_e/(rho_w-rho_e)*height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o)**3+rho_w*x_u**3/(rho_w-rho_e))**(1/3)
    return x_m

#Berechnung des Drucks

def pressure(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o):
    p_m = p_u -g_p*rho_w*(height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o)-x_u)
    return p_m

#Eigentliche Ergebnisse

def result_func(a_p, R_stern, T_stern, M_stern, R_p, A_B, e, k_2, Q, H2O):
    rho_w = 990 #mittlere Dichte von Wasser in kg/m^3
    rho_e = 910 #mittlere Dichte von Eis in kg/m^3
    
    M_H2O = rho_w*4/3*math.pi*((R_p+H2O)**3-R_p**3) #Gesamtmasse des Wassers auf dem Planeten
    M_planet = 5513*4/3*math.pi*R_p**3 #Masse des Planeten ausgerechnet über die dichte, erddichte als referenzwert
    
    #zwischenschritte
    n = math.sqrt(G*M_stern/(a_p**3))
    E_tidal = 21/2*k_2/Q*G*M_stern**2*R_p**5*n*e**2/a_p**6
    h_s = sigma*R_stern**2*T_stern**4/a_p**2
    E_stern = h_s*math.pi*R_p**2*(1-A_B)
    T_o = (E_stern/(4*math.pi*sigma*R_p**2))**(1/4)
    T_u = ((E_tidal)/(4*math.pi*sigma*R_p**2))**(1/4)

    #Naturkonstanten
    mu_w = 0.018015 #molare Masse von Wasser in kg/mol
    R = 8.314 #ideale Gaskonstante J/(mol*K)

    p_u = rho_w*R*T_u/mu_w #Druck zwischen Ozean und Gestein über ideales gas
    p_o = rho_e*R*T_o/mu_w
    
    print(p_u, p_o)
    print(p_o-p_u)
    '''
    #Temperatur in K und Dampfdruck in pa
    T = []
    p = []
    p_u = 0
    p_o = 0

    with open('Dampfdruck.csv') as tabelle_csv:
        tabelle = csv.reader(tabelle_csv, delimiter=',')
        
        zeilennummer = 0
        
        for row in tabelle:
            if zeilennummer == 0:
                zeilennummer += 1
            else: 
                T.append(float(row[0])+273.15) #Umrechnung in Celsius
                p.append(float(row[1])*100) #Umrechnung in pascal
                zeilennummer += 1

    j = 0
    k = 0

    if T_u-273.15 < 373.217:

        while T[j] < T_u:
            p_u = p[j]
            j = j+1
        
        while T[k] < T_o:
            p_o = p[k]
            k = k+1
            
    print('unterer Druck', p_u)
    print('oberer Druck', p_o)
    
    if p_u == 0:
        return 0, 0, 0, T_u
'''    
    x_u = R_p
    g_p = G*M_planet/(R_p**2) #Gravitation des Planeten
    
    x_o = height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o)
    x_m = height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o)
    
    print ('dichte',-g_p*(p_o-p_u)/(x_o-x_u))
    print (x_u-(p_o-p_u)/(g_p*990)) #alternatividee für x_o mit der Annahme dass rho_ges = rho_w, liefert etwas zu dicke eisschicht

    return x_o-x_u, x_o-x_m, x_m-x_u, T_u
        
'''
#p_m = g_p*rho_e/3*(height_of_ice(rho_e, rho_w, M_H2O, x_u)**3-height_of_water(rho_e, rho_w, M_H2O, x_u)**3)/height_of_water(rho_e, rho_w, M_H2O, x_u)**2

print(height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u))
print(height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u))
print(pressure(rho_e, rho_w, M_H2O, x_u, g_p, p_u))
print(p_u)
print(3*M_H2O/(math.pi*4*(height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u)**3-x_u**3)))
print((rho_e*(height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u)-height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u))+rho_w*(height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u)-x_u))/(height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u)-x_u))
#print(rho_w*height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u))
'''