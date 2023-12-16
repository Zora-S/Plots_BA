import math
import cmath
import csv
from scipy.constants import G
from scipy.constants import Stefan_Boltzmann as sigma

#Naturkonstanten

rho_w = 990 #mittlere Dichte von Wasser in kg/m^3
rho_e = 910 #mittlere Dichte von Eis in kg/m^3
mu_w = 0.018015 #molare Masse von Wasser in kg/mol

#angenommene Parameter

R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse
R_planet = 6e6
M_planet = 5513*4/3*math.pi*R_planet**3 #Masse des Planeten ausgerechnet über die dichte, erddichte als referenzwert 
A_B = 0.9 #0.1 bis 0.9
e = 0.5 #0.1 bis 0.9

k_2 = 0.3 #0.01 bis 1, 0.3 bei der erde
Q = 100 #1 bis 10**6

M_H2O = rho_w*4/3*math.pi*((R_planet+4000)**3-R_planet**3) #Gesamtmasse des Wassers auf dem Planeten

#zwischenschritte

a_p = 0.5 * 14959787000 #AU #math.sqrt(1-A_B)*(T_stern/273.15)**2*R_stern/2 #hier ist Bedingung 1 enthalten, mindestabstand für äußere eisschicht
n = math.sqrt(G*M_stern/(a_p**3))
E_tidal = 21/2*k_2/Q*G*M_stern**2*R_planet**5*n*e**2/a_p**6
h_s = sigma*R_stern**2*T_stern**4/a_p**2
E_stern = h_s*math.pi*R_planet**2*(1-A_B)
T_o = (E_stern/(4*math.pi*sigma*R_planet**2))**(1/4)
T_u = ((E_tidal)/(4*math.pi*sigma*R_planet**2))**(1/4)

x_u = R_planet
g_p = G*M_planet/(R_planet**2) #Gravitation des Planeten
p_u = 0 #Druck zwischen Ozean und Gestein, wird hier nur initialisiert

#Temperatur in K und Dampfdruck in pa
T = []
p = []

with open('Dampfdruck.csv') as tabelle_csv:
    tabelle = csv.reader(tabelle_csv, delimiter=',')
    
    zeilennummer = 0
    
    for row in tabelle:
        if zeilennummer == 0:
            zeilennummer += 1
        else: 
            T.append(float(row[0])+273.15)
            p.append(float(row[1])*100)
            zeilennummer += 1

j = 0

if T_u-273.15 < 373.217:

    while T[j] < T_u:
        p_u = p[j]
        j = j+1
        
else:
    print('Zu warm')
    
#Abfrage, ob der Gewichtsdruck der Gesamtwassermasse groß genug ist
if p_u < (M_H2O*g_p/(4*math.pi*R_planet**2)):
    print ('genug Wasser')
    
else:
    print ('zu wenig Wasser')

#Hilfsgrößen

a = 4/3*math.pi*rho_e
b = 4/3*math.pi*(rho_w-rho_e)
c = -x_u**3 *4/3*math.pi*rho_w - M_H2O
d = -1-rho_e/(3*rho_w)
e = x_u-p_u/(g_p*rho_w) 
f = 0
g = 0
h = rho_e/(3*rho_w)

B = e/(d-h*b/a)
C = 0
D = (-h*c/a)/(d-h*b/a)

#Bestimmung von x_i mit cardanischer formel

p = C - B**2 /3
q = 2* B**3 /27 - B*C /3 +D

Dis = p**3 /27 + q**2 /4

print ('Dis', Dis)

if Dis > 0:
    w1 = -q/2+math.sqrt(Dis)
    w2 = -q/2-math.sqrt(Dis)
    z = math.pow(abs(w1), 1/3) * (-1 if w1 < 0 else 1) + math.pow(abs(w2), 1/3) * (-1 if w2 < 0 else 1)
    
    x_m = z-B/3

    if x_m < x_u:
        print ('zu wenig wasser')

    #Bestimmung von x_o

    x_o = (-b/a*(x_m)**3-c/a)**(1/3)
    p_m = g_p*rho_e*(x_o**3-x_m**3)/(3*x_m**2)

    print ('Eisschicht:', x_o-x_m)
    print ('Wasserschicht:', x_m-x_u)
    print ('Beide Schichten', x_o-x_u)
    
elif Dis <0 :
    w1 = -q/2+cmath.sqrt(Dis)
    w2 = -q/2-cmath.sqrt(Dis)
    
    cr1_w1 = cmath.exp(cmath.log(w1) / 3)
    cr2_w1 = cr1_w1 * cmath.exp(2j * cmath.pi / 3)
    cr3_w1 = cr1_w1 * cmath.exp(4j * cmath.pi / 3)
    
    cr1_w2 = cmath.exp(cmath.log(w1) / 3)
    cr2_w2 = cr1_w2 * cmath.exp(2j * cmath.pi / 3)
    cr3_w2 = cr1_w2 * cmath.exp(4j * cmath.pi / 3)
    
    z = [cr1_w1.real + cr1_w2.real, cr2_w1.real + cr2_w2.real, cr3_w1.real + cr3_w2.real]

    x_m = []
    
    for i in range(3):
        x_m.append(z[i-1]-B/3)

    print ('x_m:', x_m)

#Bestimmung von x_o

    x_o = [] 
    
    for i in range(3):
        x_o.append((-b/a*(x_m[i-1])**3-c/a)**(1/3))
        
    print(a, b, c)
    print('x_o:', x_o)
        
    for i in range(3):
        if x_m[i-1] > x_u:
            print (x_o[i-1], x_m[i-1])
            print ('Eisschicht:', x_o[i-1]-x_m[i-1])
            print ('Wasserschicht:', x_m[i-1]-x_u)
            print ('Beide Schichten:', x_o[i-1]-x_u)

            
print ('Temperatur in °C:', T_u-273.15)
print ('Temperatur in K:', T_u)
print ('große Halbachse:', a_p)
print ('unterer Druck:', p_u)
print (M_H2O*g_p/(4*math.pi*R_planet**2))
print ('mittlerer Druck:', p_m)