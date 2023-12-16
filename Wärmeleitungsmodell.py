from scipy.optimize import fsolve
import numpy as np
import math
from scipy.constants import G
from scipy.constants import Stefan_Boltzmann as sigma

R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse
R_p = 6e6
A_B = 0.9 #0.1 bis 0.9
e = 0.05 #0.1 bis 0.9
dW = 4000 #Gesamte Schichtdicke

x_u = R_p
x_o = R_p + dW

k_2 = 0.3 #0.01 bis 1, 0.3 bei der erde
Q = 100 #1 bis 10**6

#zwischenschritte

a_p = 0.8 * 14959787000 #AU #math.sqrt(1-A_B)*(T_stern/273.15)**2*R_stern/2 #hier ist Bedingung 1 enthalten, mindestabstand für äußere eisschicht


rho_w = 997 #mittlere Dichte von Wasser in kg/m^3
rho_e = 910 #mittlere Dichte von Eis in kg/m^3
l_W = 0.5562 #Wärmeleitfähigkeit des Wassers
l_E = 2.33 #Wärmeleitfähigkeit des Eises

M_planet = 5513*4/3*math.pi*R_p**3 #Masse des Planeten ausgerechnet über die dichte, erddichte als referenzwert
g_p = G*M_planet/(R_p**2) #Gravitation des Planeten
    
    #zwischenschritte
n = math.sqrt(G*M_stern/(a_p**3))
E_tidal = 21/2*k_2/Q*G*M_stern**2*R_p**5*n*e**2/a_p**6
h_s = sigma*R_stern**2*T_stern**4/a_p**2
E_stern = h_s*math.pi*R_p**2*(1-A_B)
T_o = ((E_stern+E_tidal)/(4*math.pi*sigma*R_p**2))**(1/4) #Berechnung mit oder ohne leistung der Gezeitenheizung?
T_m = 273.15 #mittlere Temperatur wird auf den Gefirerpunkt gesetzt
A = 4*math.pi*x_u**2 #Oberfläche des festen Kerns in Quadratmetern

x_m = x_o - (l_E*A) * (T_m-T_o)/(E_tidal)

T_u = T_m + E_tidal/(l_W*A) * (x_m-x_u)

print("große Halbachse in AU", a_p/14959787000)
print("untere Leistung pro Quadratmeter", E_tidal/A)
print("obere Leistung pro Quadratmeter", E_stern/A)
print("obere Temperatur in Kelvin", T_o)
print("untere Temperatur in Kelvin", T_u)
print("Dicke der Wasserschicht in Metern", x_m-x_u)
print("Dicke der Eisschicht in Metern", x_o-x_m)