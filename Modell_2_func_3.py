from scipy.optimize import fsolve
import numpy as np
import math
import csv
from scipy.constants import G
from scipy.constants import Stefan_Boltzmann as sigma

# Gleichung: 0 = a*x**3 + b + (c + d*(x-y)/(x**3-y**3))*(e + f*x**3)**(2/3)

R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse
R_p = 6e6
A_B = 0.9 #0.1 bis 0.9
e = 0.5 #0.1 bis 0.9
#M_planet = 5513*4/3*math.pi*R_p**3

x_u = R_p

k_2 = 0.3 #0.01 bis 1, 0.3 bei der erde
Q = 100 #1 bis 10**6

H2O = 10000

#zwischenschritte

a_p = 0.5 * 14959787000 #AU #math.sqrt(1-A_B)*(T_stern/273.15)**2*R_stern/2 #hier ist Bedingung 1 enthalten, mindestabstand für äußere eisschicht


rho_w = 997 #mittlere Dichte von Wasser in kg/m^3
rho_e = 910 #mittlere Dichte von Eis in kg/m^3
l_W = 0.5562 #wärmeleitfähigkeit des wassers
    
M_H2O = rho_w*4/3*math.pi*((R_p+H2O)**3-R_p**3) #Gesamtmasse des Wassers auf dem Planeten
M_planet = 5513*4/3*math.pi*R_p**3 #Masse des Planeten ausgerechnet über die dichte, erddichte als referenzwert
g_p = G*M_planet/(R_p**2) #Gravitation des Planeten
    
    #zwischenschritte
n = math.sqrt(G*M_stern/(a_p**3))
E_tidal = 21/2*k_2/Q*G*M_stern**2*R_p**5*n*e**2/a_p**6
h_s = sigma*R_stern**2*T_stern**4/a_p**2
E_stern = h_s*math.pi*R_p**2*(1-A_B)
T_o = (E_stern/(4*math.pi*sigma*R_p**2))**(1/4)
T_u = 273.15+50 #((E_tidal)/(4*math.pi*sigma*R_p**2))**(1/4)

x_m = l_W*4*math.pi*x_u**2*(T_u-273.15)/(E_tidal) #nur wasserschicht, ohne planetenkern
x = 500
T_u2 = E_tidal*x/(l_W*4*math.pi*x_u**2) - 273.15

print(T_u)
print(T_u2)
print(x_m)

'''
    #Naturkonstanten
mu_w = 0.018015 #molare Masse von Wasser in kg/mol
R = 8.314 #ideale Gaskonstante J/(mol*K)

p_u = rho_w*g_p*H2O #Druck zwischen Ozean und Gestein über ideales gas und gewichtsdruck
p_o = 100 #rho_e*R*T_o/mu_w

print('obere Temperatur/K', T_o)
print('untere Temperatur/K', T_u)
print('Druckdifferenz', p_o)
print('p_u in bar', p_u/(10**5))
print('druck ideales gas',rho_w*R*T_u/mu_w/(10**5) )

x_o = x_u + H2O

x_m = - rho_e/(rho_w-rho_e)*x_o - p_o/(g_p*(rho_w-rho_e)) + rho_w/(rho_w-rho_e)*x_u + p_u/(g_p*(rho_w-rho_e))

print(x_m)
'''
'''
def height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o):
    # Parameterwerte
    a = (rho_w-rho_e)*g_p
    b = g_p*rho_e
    c = (rho_e-rho_w)/rho_e
    d = 3*M_H2O/(4*math.pi*rho_e)+(rho_w/rho_e *x_u**3)
    e = p_o-p_u-g_p*rho_w*x_u

    # Gleichung definieren
    equation = lambda x: a*x+b*(c*x**3+d)**(1/3)+e
    diff = lambda x: a+((b*c*x**2)/(c*x**3+d)**(2/3))

    # Schätzwert für die Lösung
    initial_guess = x_u+6000
    
    x_m2 = 0
    
    while (x_m-x_m2) > 1:
        x_m2 = x_m - equation(x_m)/diff(x_m)
        x_m = x_m2
        print(x_m)
        

    # Anwenden des Newton-Verfahrens
    solution = fsolve(equation, initial_guess)
    
    return solution

    # Ergebnis ausgeben
    #print("Die Lösung ist x =", solution[0])

print(height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o))
print(height_of_water(rho_e, rho_w, M_H2O, x_u, g_p, p_u, p_o)-x_u) '''