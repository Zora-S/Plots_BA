import math
import csv
from scipy.optimize import fsolve
import numpy as np
from scipy.constants import G
from scipy.constants import Stefan_Boltzmann as sigma
from Modell_2_func_2 import height_of_ice

a = 1
b = 0
c = 1
d = 1

    # Koeffizienten in einer Liste, um sie an np.roots zu übergeben
coefficients = [a, b, c, d]

    # Berechnen der Wurzeln
roots = np.roots(coefficients)
print(roots)

print(1/4*3)

#print(height_of_ice(0, 0, math.pi, 1, 2, 2))
#height_of_ice(rho_e, rho_w, M_H2O, x_u, g_p, p_u)
