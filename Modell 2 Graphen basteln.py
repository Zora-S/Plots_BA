import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from Modell_2_as_function import mindestabstand
from Modell_2_as_function import maximalabstand
from Modell_2_as_function import wasserschicht_planet



a_min = [] #große Halbachse mindestens
a_max = [] # große Halbachse maximal


#graphen für sternabhängigkeit
R_stern = [6,3.5,1.8,1.5,1,0.8,0.6,0.1] 
T_stern = [30e3,10e3,7.5e3,6e3,5.3e3,3.9e3,2.55e3,2.25e3] # in Kelvin
M_stern = [18,5,1.9,1.4,1,0.8,0.6,0.1] 
R_planet = 10e6
M_planet = 5513*4/3*math.pi*R_planet**3 #Masse des Planeten ausgerechnet über die dichte, erddichte als referenzwert 
A_B = 0.9 #0.1 bis 0.9
e = 0.5 #0.1 bis 0.9

k_2 = 0.3 #0.01 bis 1, 0.3 bei der erde
Q = 100 #1 bis 10**6

for i in range(len(R_stern)):
    R_stern[i-1] = R_stern[i-1]*6.96342*10**8 # *6.96342*10**8 weil sonnenradius
    M_stern[i-1] = M_stern[i-1]*1.989*10**30 # *1.989*10**30 weil sonnenmasse
    
for i in range(len(R_stern)):
    a_min.append(mindestabstand(R_stern[i], T_stern[i], A_B))
    a_max.append(maximalabstand(M_stern[i],273.15, R_planet, e, k_2, Q))

x1 = [i * 1/149597870000 for i in a_min]
x2 = [i * 1/149597870000 for i in a_max]
y = T_stern

f=lambda x,m,b: m*x+b 

fig, ax = plt.subplots()
ax.plot(x1,y,'rx',label="Mindestabstand, bestimmt durch Heizleistung des Sterns")
ax.plot(x2,y,'bx',label="Maximalabstand, bestimmt durch Gezeitenheizung")
ax.legend()
ax.set_xlim(0,1.1)
ax.set_ylim(1800,7700)
ax.grid(True)
ax.set_xlabel("große Halbachse in AU")
ax.set_ylabel("Temperatur des Sterns in K")

plt.show()

plt.savefig('Abschätzung ohne O,B-Stern.pdf')

'''

# Graphen für Planetenabhängigkeit
a_s = [] #Abstand durch Strahlung
R_p = [2e6, 2.5e6, 3e6, 3.5e6, 4e6, 4.5e6, 5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6, 8e6, 8.5e6, 9e6, 9.5e6, 10e6, 10.5e6, 11e6, 11.5e6, 12e6]

R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse 
A_B = 0.9
e = 0.5
k_2 = 0.3
Q = 100


for i in range(0,len(R_p)):
    
    if mindestabstand(R_stern, T_stern, A_B) < maximalabstand(M_stern, 423.15, R_p[i], e, k_2, Q):
        a_min.append(maximalabstand(M_stern, 423.15, R_p[i], e, k_2, Q))
        
    else:
        a_min.append(-1e9)
        
    if mindestabstand(R_stern, T_stern, A_B) > maximalabstand(M_stern, 423.15, R_p[i], e, k_2, Q):
        a_s.append(mindestabstand(R_stern, T_stern, A_B))
        
    else:
        a_s.append(-1e9)
        
    a_max.append(maximalabstand(M_stern, 223.15, R_p[i], e, k_2, Q))
    
print (len(a_min), len(a_max), len(R_p))
    
x1 = [i * 1/149597870000 for i in a_min]
x2 = [i * 1/149597870000 for i in a_max]
x3 = [i * 1/149597870000 for i in a_s]
y = R_p

f=lambda x,m,b: m*x+b 

rcParams["font.size"]=16
fig, ax = plt.subplots()
ax.plot(x1,y,'rx',label="$T_u$ = 150°C")
ax.plot(x2,y,'bx',label="$T_u$ = -50°C")
ax.plot(x3,y,'gx',label="Mindestabstand durch \n Strahlung des Sterns")
ax.legend(loc='lower right')
ax.set_xlim(0.02, 0.1)
ax.set_ylim(1.5e6,12.5e6)
ax.grid(True)
ax.set_xlabel("große Halbachse in AU")
ax.set_ylabel("Radius des Planeten in m")

plt.show()
fig.savefig('Abhängigkeit_Planetenradius_mit_e=0,5.pdf')



a_s = [] #Abstand durch Strahlung
R_p = 10e6
R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse 
A_B = 0.9
e = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
k_2 = 0.3
Q = 100


for i in range(0,len(e)):
    
    if mindestabstand(R_stern, T_stern, A_B) < maximalabstand(M_stern, 423.15, R_p, e[i], k_2, Q):
        a_min.append(maximalabstand(M_stern, 423.15, R_p, e[i], k_2, Q))
        
    else:
        a_min.append(-1e9)
        
    if mindestabstand(R_stern, T_stern, A_B) > maximalabstand(M_stern, 423.15, R_p, e[i], k_2, Q):
        a_s.append(mindestabstand(R_stern, T_stern, A_B))
        
    else:
        a_s.append(-1e9)
        
    a_max.append(maximalabstand(M_stern, 223.15, R_p, e[i], k_2, Q))
    
print (len(a_min), len(a_max), len(e))

#Umrechnung in AU    
x1 = [i * 1/149597870000 for i in a_min]
x2 = [i * 1/149597870000 for i in a_max]
x3 = [i * 1/149597870000 for i in a_s]
y = e

f=lambda x,m,b: m*x+b 

#rcParams["figure.figsize"]=6.3 ,10
rcParams["font.size"]=16
fig, ax = plt.subplots()
ax.plot(x1,y,'rx',label="$T_u$ = 150°C")
ax.plot(x2,y,'bx',label="$T_u$ = -50°C")
ax.plot(x3,y,'gx',label="Mindestabstand durch \n Strahlung des Sterns")
ax.legend(loc='lower right')
yticks = np.arange(0.1, 1, 0.1)  # Start, Stop, Schrittweite
plt.yticks(yticks) #Setzen der Beschriftung
ax.set_xlim(0.03, 0.1)
#ax.set_ylim(1.5e6,12.5e6)
ax.grid(True)
ax.set_xlabel("große Halbachse in AU")
ax.set_ylabel("Exzentrizität der Bahn")

plt.show()
#fig.savefig('Abhängigkeit_Exzentrizität.pdf')



#Wassermasse plot
R_p = 8e6

a_p = 0.05 * 149597870000
R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse 
A_B = 0.9
e = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
k_2 = 0.3
Q = 100
H2O = 4000

beide = []
eis = []
wasser = []
temp = []

for i in range(len(e)):
    result = wasserschicht_planet(a_p, R_stern, T_stern, M_stern, R_p, A_B, e[i], k_2, Q, H2O)
    #print (result)
    beide.append(result[0]) #dicke beider Schichten in m
    eis.append(result[1]) #dicke der Eisschicht in m
    wasser.append(result[2]) #dicke der Wasserschicht in m
    temp.append(result[3])
    
#print(beide)
    
x2 = []
x3 = []
    
for i in range(len(e)):
    x2.append(eis[i]/(beide[i]+1)*100)
    x3.append(wasser[i]/(beide[i]+1)*100)

y = temp
#y2 = temp

f=lambda x,m,b: m*x+b 

rcParams["font.size"]=16
fig, ax = plt.subplots()
#ax.plot(x1,y,'rx',label="Dicke beider Schichten")
ax.plot(x2,y,'bx',label="Dicke der Eisschicht")
ax.plot(x3,y,'gx',label="Dicke der Wasserschicht")
ax.legend(loc='lower left')
ax.grid(True)
ax.set_xlim(0, 100)
ax.set_ylim(180, 620)
ax.set_xlabel("Schichtdicken in Prozent der Gesamtdicke")
ax.set_ylabel("Temperatur in K")

ax2 = ax.twinx()

# Plotten Sie die Daten auf der zweiten Y-Achse
ax2.plot(x1,y2,'rx',label="Dicke beider Schichten")
ax2.plot(x2,y2,'bx',label="Dicke der Eisschicht")
ax2.plot(x3,y2,'gx',label="Dicke der Wasserschicht")
#ax2.grid(True)
#ax2.set_ylim(260,630)
ax2.set_ylabel("Temperatur in K", fontsize=20)

# Fügen Sie eine Legende für beide Y-Achsen hinzu
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

plt.show()
fig.savefig('Abhängigkeit_Schichtdicke_von_Temperatur_und_Exzentrizität.pdf')
'''