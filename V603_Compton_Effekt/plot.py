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

E_Ka = 8.038*10**(3)*scipy.constants.e
E_Kb = 8.905*10**(3)*scipy.constants.e
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

print("Spektum von Kupfer")

theta_Cu,N = np.genfromtxt("EmissionCu.dat", unpack = True)
plt.plot(theta_Cu,N,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
N_loc = scipy.signal.find_peaks(N,height=1000)
N_peak = np.array([1599,5050])
theta_peak = np.array([theta_Cu[122],theta_Cu[145]])
plt.plot(theta_peak,N_peak,'k|',label = 'charakteristische Peaks')
plt.plot(theta_Cu[1:119],N[1:119],'b--',label = 'Bremsberg')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/CU_Spektrum.pdf",bbox_inches='tight')
plt.close()

print("Berechnung der Energien der Peaks")

theta_peak= theta_peak*((360)**(-1))*(2*np.pi)
lambda_peak = ([(2*d_LiF)*np.sin(theta_peak[0]), (2*d_LiF)*np.sin(theta_peak[1])])
E_peak = ([(h*c)/lambda_peak[0], (h*c)/lambda_peak[1]])
E_peak= ([E_peak[0]*6.242*10**(18),E_peak[1]*6.242*10**(18)])

print("Transmission")

theta_ComptonOhne,N_0 = np.genfromtxt("ComptonOhne.txt", unpack = True)
theta_ComptonAl,N_Al = np.genfromtxt("ComptonAl.txt", unpack = True)
N_0 = unp.uarray(N_0,np.sqrt(N_0))
N_Al = unp.uarray(N_Al,np.sqrt(N_Al))

print("Totzeitkorrektur")

I_0 = (N_0)/(1-(90*10**(-6)* N_0))
I_Al = (N_Al)/(1-(90*10**(-6)* N_Al))
T = I_Al/I_0
theta_ComptonOhne = theta_ComptonOhne*((360)**(-1))*(2*np.pi)
lambda_Compton = 2*d_LiF*np.sin(theta_ComptonOhne)


print("Plot und Ausgleichskurve")
def f(x,m,b):
    return m * x + b

params,cov = curve_fit(f,lambda_Compton/10**(-12),noms(T))
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)
plt.plot(lambda_Compton,noms(T),'x',label = 'Messdaten')
plt.plot(lambda_Compton,f(lambda_Compton/10**(-12),*params),'r-',label = 'Ausgleichsgerade')
plt.ylabel(f"Transmission")
plt.xlabel(f"Wellenlänge in Meter")
plt.legend()
plt.savefig("plots/Transmission.pdf",bbox_inches='tight')
plt.close()

print("Compton Wellenlänge")

I = unp.uarray(2731,np.sqrt(2731))
I_1 = unp.uarray(1180,np.sqrt(1180))
I_2 = unp.uarray(1024,np.sqrt(1024))

T_1 = I_1/I
T_2 = I_2/I

lambda_1 = (T_1-unparams[1])/unparams[0]
lambda_2 = (T_2-unparams[1])/unparams[0]
lambda_C_berechnet = lambda_2 - lambda_1



print(f"""
Vorbereitung:
lambda A ist:{lambda_Ka}
lambda B ist:{lambda_Kb}
alpha A ist:{alpha_Ka}
alpha B ist:{alpha_Kb}
Compton Wellenlänge ist:{lambda_C}
und in der Literatur: {lambda_Clit}

Kupferspektrum:
Die Peaks sind bei den Winkeln: {theta_peak}
mit den Impulsen : {N_peak}
Die Wellenlängen an den Peaks ist: {lambda_peak}
Die Energie an den Peaks ist: {E_peak}

Transmission:
Parameter der Ausgleichsgeraden:
m ist {unparams[0]} und b ist {unparams[1]}

Compton Wellenlänge:
I1 ist :{I_1}
I2 ist :{I_2}
I ist :{I}
T1 ist : {T_1}
T2 ist : {T_2}
lambda 1 [pm]:{lambda_1}
lambda 2 [pm]:{lambda_2}
Comptonwellenlänge in Pikometer : {lambda_C_berechnet}
""")










