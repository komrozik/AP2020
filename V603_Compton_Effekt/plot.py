import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
#np.seterr(divide='ignore', invalid='ignore')
import sympy as sym
from scipy import integrate
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.optimize import curve_fit
import scipy.constants

#Vorbereitung

E_Ka = 8.038*10**(3)
E_Kb = 8.905*10**(3)
h = scipy.constants.h
c = scipy.constants.c
m_e = scipy.constants.m_e
d_LiF = 201.4*10**(-12)
lambda_Ka = (h*c)/E_Ka
lambda_Kb = (h*c)/E_Kb
lambda_C = h/(c*m_e)
lambda_Clit = scipy.constants.physical_constants['Compton wavelength']
alpha_Ka = np.arcsin((lambda_Ka)/(2*d_LiF))
alpha_Kb = np.arcsin((lambda_Kb)/(2*d_LiF))
print(f'lambda A ist:{lambda_Ka}')
print(f'lambda B ist:{lambda_Kb}')
print(f'alpha A ist:{alpha_Ka}')
print(f'alpha B ist:{alpha_Kb}')
print(f'Compton Wellenlänge ist:{lambda_C} und in der Literatur{lambda_Clit}: ')









