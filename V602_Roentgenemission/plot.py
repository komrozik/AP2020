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
R_y = scipy.constants.value("Rydberg constant times hc in J")
d_LiF = 201.4*10**(-12)
Z = np.array([30,32,35,37,38,40])
E_K = np.array([9.65,11.10,13.47,15.20,16.10,17.99])
E_K = E_K*scipy.constants.e*10**(3)
lambda_K = (h*c)/E_K
theta_K = np.arcsin((1*lambda_K)/(2*d_LiF))
theta_K = (theta_K*180)/(np.pi)
sigma_K = Z-np.sqrt((E_K)/(R_y)-((alpha**2)*(Z**4))/(4))
lambda_Ka = (h*c)/E_C_Ka
theta_Ka = np.arcsin((1*lambda_Ka)/(2*d_LiF))
theta_Ka = (theta_Ka*180)/(np.pi)
lambda_Kb = (h*c)/E_C_Kb
theta_Kb = np.arcsin((1*lambda_Kb)/(2*d_LiF))
theta_Kb = (theta_Kb*180)/(np.pi)


print("Braggbedingung")

theta_Bragg,N_Bragg = np.genfromtxt("data/Bragg.dat", unpack = True)
plt.plot(theta_Bragg,N_Bragg,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
N_loc = scipy.signal.find_peaks(N_Bragg,height=100)
#   Peak bei 22
N_peak = N_Bragg[N_loc[0]]
theta_peak = theta_Bragg[N_loc[0]]
plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Bragg.pdf",bbox_inches='tight')
plt.close()
#Vergleich Bragg winkel:
theta_abs = np.sqrt((28-theta_peak)**2)
theta_rel = np.sqrt(((theta_abs)/28)**2)

print("Emissionsspektrum")

theta_Cu,N_Cu = np.genfromtxt("data/Emissionsspektrum.dat", unpack = True)
plt.plot(theta_Cu,N_Cu,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
N_loc = scipy.signal.find_peaks(N_Cu,height=1000)
#   Peak bei 22
peak_loc = N_loc[0]
N_peak_Cu = N_Cu[N_loc[0]]
theta_peak_Cu = theta_Cu[N_loc[0]]
plt.plot(theta_peak_Cu,N_peak_Cu,'k|',label = 'Peaks')
plt.xlabel(f"Winkel / Grad")
plt.plot(theta_Cu[0:peak_loc[0]-2],N_Cu[0:peak_loc[0]-2],'b--',label = 'Bremsberg')
plt.plot(theta_Cu[peak_loc[0]+5:peak_loc[1]-4],N_Cu[peak_loc[0]+5:peak_loc[1]-4],'b--')
plt.plot(theta_Cu[peak_loc[1]+6:180],N_Cu[peak_loc[1]+6:180],'b--')
plt.legend()
plt.savefig("plots/Cu_Emission.pdf",bbox_inches='tight')
plt.close()

print("Absorptionskanten")

#Brom
theta_Brom,N_Brom = np.genfromtxt("data/Brom.dat", unpack = True)
plt.plot(theta_Brom,N_Brom,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Brom,height=100)
# #   Peak bei 22
# N_peak = N_Brom[N_loc[0]]
# theta_peak = theta_Brom[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Brom.pdf",bbox_inches='tight')
plt.close()

#Gallium
theta_Gallium,N_Gallium = np.genfromtxt("data/Gallium.dat", unpack = True)
plt.plot(theta_Gallium,N_Gallium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Gallium,height=100)
# #   Peak bei 22
# N_peak = N_Gallium[N_loc[0]]
# theta_peak = theta_Gallium[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Gallium.pdf",bbox_inches='tight')
plt.close()

#Rubidium
theta_Rubidium,N_Rubidium = np.genfromtxt("data/Rubidium.dat", unpack = True)
plt.plot(theta_Rubidium,N_Rubidium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Rubidium,height=100)
# #   Peak bei 22
# N_peak = N_Rubidium[N_loc[0]]
# theta_peak = theta_Rubidium[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Rubidium.pdf",bbox_inches='tight')
plt.close()

#Strontium
theta_Strontium,N_Strontium = np.genfromtxt("data/Strontium.dat", unpack = True)
plt.plot(theta_Strontium,N_Strontium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Strontium,height=100)
# #   Peak bei 22
# N_peak = N_Strontium[N_loc[0]]
# theta_peak = theta_Strontium[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Strontium.pdf",bbox_inches='tight')
plt.close()

#Zink
theta_Zink,N_Zink = np.genfromtxt("data/Zink.dat", unpack = True)
plt.plot(theta_Zink,N_Zink,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Zink,height=100)
# #   Peak bei 22
# N_peak = N_Zink[N_loc[0]]
# theta_peak = theta_Zink[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Zink.pdf",bbox_inches='tight')
plt.close()

#Zirkonium
theta_Zirkonium,N_Zirkonium = np.genfromtxt("data/Zirkonium.dat", unpack = True)
plt.plot(theta_Zirkonium,N_Zirkonium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
# N_loc = scipy.signal.find_peaks(N_Zirkonium,height=100)
# #   Peak bei 22
# N_peak = N_Zirkonium[N_loc[0]]
# theta_peak = theta_Zirkonium[N_loc[0]]
# plt.plot(theta_peak,N_peak,'k|',label = 'Maximum')
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Zirkonium.pdf",bbox_inches='tight')
plt.close()


print(f"""
VORBEREITUNG:
Rydbergenergie: {R_y} Joule
Feinstruktur: {alpha}
Energien: {E_K} Joule
Vorbereitung:
Die Bragg Winkel bei den K-Kanten, der Beugungsordnung 1 und mit LiF Kristall sind:
{theta_K} in °
Die Abschirmkonstante ist jeweils:
{sigma_K}
Für Kupfer ist:
K alpha: {E_C_Ka}
K beta: {E_C_Kb}
Mit den Bragg Winkeln:
Alpha: {theta_Ka}
Beta: {theta_Kb}
----------------------
BRAGGBEDINGUNG:
Das Maximum der Kurve weicht um {theta_abs}° vom Sollwert ab.
Das ist ein relativer Fehler von {theta_rel}. Ungefähr {theta_rel*100} %.

""")
