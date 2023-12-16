from Modell_2_func_2 import result_func
import math
from scipy.constants import G

#angenommene Parameter

R_stern = 0.6 *6.96342*10**8 # *6.96342*10**8 weil sonnenradius
T_stern = 2500 # in Kelvin
M_stern = 0.6 *1.989*10**30 # *1.989*10**30 weil sonnenmasse
R_p = 6e6
A_B = 0.9 #0.1 bis 0.9
e = 0.5 #0.1 bis 0.9
#M_planet = 5513*4/3*math.pi*R_p**3

k_2 = 0.3 #0.01 bis 1, 0.3 bei der erde
Q = 100 #1 bis 10**6

H2O = 400000

#zwischenschritte

a_p = 0.5 * 14959787000 #AU #math.sqrt(1-A_B)*(T_stern/273.15)**2*R_stern/2 #hier ist Bedingung 1 enthalten, mindestabstand für äußere eisschicht

result = result_func(a_p, R_stern, T_stern, M_stern, R_p, A_B, e, k_2, Q, H2O)

#print(990*(result[2]+R_p))
print('Beide Schichten', result[0])
print('Eisschicht', result[1])
print('Wasserschicht', result[2])
print('Temperatur', result[3])
print('Celsius', result[3]-273.15)