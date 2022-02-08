#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:13:26 2022

@author: Raoul
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.size'] = 11

delta_nu = 50 # just a first guess to make the plot look nice
n = np.arange(39, 47) # just a first guess to make the plot look nice
D = 0.75 # just a first guess to make the plot look nice
def nu(delta_nu, n, l):
    return delta_nu * (n + l/2 + 1/4) - D * l * (l + 1)

# define frequency arrays for each value of l
nu_0 = nu(delta_nu, n, 0)
nu_1 = nu(delta_nu, n, 1)
nu_2 = nu(delta_nu, n, 2)
nu_3 = nu(delta_nu, n, 3)

def lorentzian(x, x_0, gamma):
    return 1 / (np.pi * gamma * (1 + ((x - x_0) / gamma)**2))

gamma = 0.5 # just a first guess to make the plot look nice
x = np.linspace(2100, 2400, 10000) # range of frequencies to be plotted

plt.figure(figsize=(8, 5))

Amp_solar = 2.1 #ppm 
L_solar = 3.846e26 #Watts 
M_solar = 1.989e30 #kg
Teff_solar = 5777 #K
deltaT = 1250 #K

L = L_solar
M = M_solar
T_eff = Teff_solar

Tred = 8907*((L/L_solar)**-0.093)
beta = 1 - np.exp((T_eff - Tred)/deltaT)
amplitude = Amp_solar*beta*(L/L_solar)*(M_solar/M)*(Teff_solar/T_eff)**2

print(amplitude)

T0 = 436 #+- 24 K
gamma0 = 1.02 #+- 0.07 microHz

#Linewidth
linewidth = gamma0*np.exp((T_eff-Teff_solar)/T0)

#Height in FREQUENCYDOMAIN
Height = (2*amplitude**2)/(np.pi*linewidth)


# use for loops to plot all the values of n and l
for i in nu_0:
    y = 7.5 * lorentzian(x, i, gamma)
    plt.plot(x, y)
for i in nu_1:
    y = 10 * lorentzian(x, i, gamma)
    plt.plot(x, y)
for i in nu_2:
    y = 4 * lorentzian(x, i, gamma)
    plt.plot(x, y)
for i in nu_3:
    y = lorentzian(x, i, gamma)
    plt.plot(x, y)
    
plt.title('Oscillation spectrum sample plot')
plt.xlabel(r'Frequency [$\mu$Hz]')
plt.ylabel('PSD')
plt.savefig('oscillation_spectrum.pdf')
