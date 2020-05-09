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
from scipy.signal import find_peaks
import scipy.constants

print("Vorbereitung")

E_C_Ka = 8.038*10**(3)*scipy.constants.e
E_C_Kb = 8.905*10**(3)*scipy.constants.e
h = scipy.constants.h
c = scipy.constants.c
m_e = scipy.constants.m_e
alpha = scipy.constants.value("fine-structure constant")
d_LiF = 201.4*10**(-12)
Z = np.array([30,32,35,37,38,40])
E_K = np.array([9.65,11.10,13.47,15.20,16.10,17.99])
E_K = E_K*scipy.constants.e*10**(3)
lambda_K = (h*c)/E_K
theta_K = np.arcsin((1*lambda_K)/(2*d_LiF))
theta_K = (theta_K*180)/(np.pi)


print(f"""
Vorbereitung:
Die Bragg Winkel bei den K-Kanten, der Beugungsordnung 1 und mit LiF Kristall sind:
{theta_K} in °
{alpha}
""")
