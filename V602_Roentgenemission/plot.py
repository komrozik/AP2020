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
from scipy.interpolate import UnivariateSpline

print("Vorbereitung")

E_C_Ka = 8.038*10**(3)*scipy.constants.e
E_C_Kb = 8.905*10**(3)*scipy.constants.e
h = scipy.constants.h
c = scipy.constants.c
m_e = scipy.constants.m_e
alpha = scipy.constants.value("fine-structure constant")
R_y = scipy.constants.value("Rydberg constant times hc in J")
d_LiF = 201.4*10**(-12)
Z = np.array([30,31,32,35,37,38,40])
E_K = np.array([9.65,10.37,11.10,13.47,15.20,16.10,17.99])
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

plt.plot(theta_Cu[peak_loc[0]-15:peak_loc[1]+15],N_Cu[peak_loc[0]-15:peak_loc[1]+15],'rx',label = 'Messdaten')
plt.plot(theta_peak_Cu,N_peak_Cu,'k|',label = 'Peaks')
spline = UnivariateSpline(theta_Cu[peak_loc[0]+5:peak_loc[1]+15], N_Cu[peak_loc[0]+5:peak_loc[1]+15]-np.max(N_Cu[peak_loc[0]-15:peak_loc[1]+15])/2, s=0)
r1a, r2a = spline.roots()
plt.axvspan(r1a, r2a, facecolor='g', alpha=0.5)
spline = UnivariateSpline(theta_Cu[peak_loc[0]-15:peak_loc[1]-5], N_Cu[peak_loc[0]-15:peak_loc[1]-5]-np.max(N_Cu[peak_loc[0]-15:peak_loc[1]-5])/2, s=0)
r1b, r2b = spline.roots() 
plt.axvspan(r1b, r2b, facecolor='g', alpha=0.5)
plt.ylabel(f"Impulse pro Sekunde")
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Cu_Peaks.pdf",bbox_inches='tight')
plt.close()

FWHMa = (h*c)/(2*d_LiF*np.sin(r1a*np.pi/180))-(h*c)/(2*d_LiF*np.sin(r2a*np.pi/180))
FWHMb = (h*c)/(2*d_LiF*np.sin(r1b*np.pi/180))-(h*c)/(2*d_LiF*np.sin(r2b*np.pi/180))
E_K_peaka= (h*c)/(2*d_LiF*np.sin(theta_peak_Cu[1]*np.pi/180))
E_K_peakb= (h*c)/(2*d_LiF*np.sin(theta_peak_Cu[0]*np.pi/180))
Aa = E_K_peaka/FWHMa
Ab = E_K_peakb/FWHMb

#Sigmas Bestimmen
E_K_lit = 8.98*10**3*scipy.constants.e
sigma_1 = 29-np.sqrt((E_K_lit)/(R_y))
sigma_2 = 29-np.sqrt( 4*(29-sigma_1)**2-(E_K_peaka)/(R_y)*4 )
sigma_3 = 29-np.sqrt( 9*(29-sigma_1)**2-(E_K_peakb)/(R_y)*9 )
print("Absorptionskanten")

#Zink
theta_Zink,N_Zink = np.genfromtxt("data/Zink.dat", unpack = True)
plt.plot(theta_Zink,N_Zink,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Zink.pdf",bbox_inches='tight')
plt.close()
I_K_Zink = N_Zink[5]+(N_Zink[10]-N_Zink[5])/2
E_K_Zink = (h*c)/(2*d_LiF*np.sin(18.7*np.pi/180))
sigma_K_Zink = 30-np.sqrt((E_K_Zink)/(R_y)-((alpha**2)*(30**4))/(4))

#Gallium
theta_Gallium,N_Gallium = np.genfromtxt("data/Gallium.dat", unpack = True)
plt.plot(theta_Gallium,N_Gallium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Gallium.pdf",bbox_inches='tight')
plt.close()
I_K_Gallium = N_Gallium[1]+(N_Gallium[6]-N_Gallium[1])/2
E_K_Gallium = (h*c)/(2*d_LiF*np.sin(17.375*np.pi/180))
sigma_K_Gallium = 31-np.sqrt((E_K_Gallium)/(R_y)-((alpha**2)*(31**4))/(4))

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
I_K_Brom = N_Brom[2]+(N_Brom[7]-N_Brom[2])/2
E_K_Brom = (h*c)/(2*d_LiF*np.sin(13.2*np.pi/180))
sigma_K_Brom = 35-np.sqrt((E_K_Brom)/(R_y)-((alpha**2)*(35**4))/(4))

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
I_K_Rubidium = N_Rubidium[3]+(N_Rubidium[9]-N_Rubidium[3])/2
E_K_Rubidium = (h*c)/(2*d_LiF*np.sin(11.8*np.pi/180))
sigma_K_Rubidium = 37-np.sqrt((E_K_Rubidium)/(R_y)-((alpha**2)*(37**4))/(4))

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
I_K_Zirkonium = N_Zirkonium[0]+(N_Zirkonium[7]-N_Zirkonium[0])/2
E_K_Zirkonium = (h*c)/(2*d_LiF*np.sin(9.95*np.pi/180))
sigma_K_Zirkonium = 40-np.sqrt((E_K_Zirkonium)/(R_y)-((alpha**2)*(40**4))/(4))

#Strontium
theta_Strontium,N_Strontium = np.genfromtxt("data/Strontium.dat", unpack = True)
plt.plot(theta_Strontium,N_Strontium,'rx',label = 'Messdaten')
plt.ylabel(f"Impulse pro Sekunde")
plt.xlabel(f"Winkel / Grad")
plt.legend()
plt.savefig("plots/Strontium.pdf",bbox_inches='tight')
plt.close()
I_K_Strontium = N_Strontium[4]+(N_Strontium[9]-N_Strontium[4])/2
E_K_Strontium = (h*c)/(2*d_LiF*np.sin(11.1*np.pi/180))
sigma_K_Strontium = 38-np.sqrt((E_K_Strontium)/(R_y)-((alpha**2)*(38**4))/(4))


def f(x,m,b):
    return m*x+b

Z = np.array([30,31,35,37,38,40])
E_K_exp = np.array([E_K_Zink,E_K_Gallium,E_K_Brom,E_K_Rubidium,E_K_Strontium,E_K_Zirkonium])
sigma_K_exp = np.array([sigma_K_Zink,sigma_K_Gallium,sigma_K_Brom,sigma_K_Rubidium,sigma_K_Strontium,sigma_K_Zirkonium])

plt.plot(Z,np.sqrt(E_K_exp),'rx',label = 'Messdaten')
params, cov = curve_fit(f,Z,np.sqrt(E_K_exp))
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)
plt.plot(Z,f(Z,*params),'b--',label = 'Ausgleichsgerade')
plt.ylabel(f"xxx")
plt.xlabel(f"xxx")
plt.legend()
plt.savefig("plots/Rydberg.pdf",bbox_inches='tight')
plt.close()
Ryd_exp = unparams[0]**2/h

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
Im Experiment:
theta: {theta_peak_Cu}

Die Sigmas:
Sigma 1: {sigma_1}
Sigma 2: {sigma_2}
Sigma 3: {sigma_3}
----------------------
BRAGGBEDINGUNG:
Das Maximum der Kurve weicht um {theta_abs}° vom Sollwert ab.
Das ist ein relativer Fehler von {theta_rel}. Ungefähr {theta_rel*100} %.

Die Auflösungen sind:
Energie a:{E_K_peaka*6.242*10**18}
Energie b:{E_K_peakb*6.242*10**18}
Delta a : {FWHMa*6.242*10**18}
Delta b: {FWHMb*6.242*10**18}
Energie 1 a:{(h*c)/(2*d_LiF*np.sin(r1a*np.pi/180))*6.242*10**18}
Energie 2 a:{(h*c)/(2*d_LiF*np.sin(r2a*np.pi/180))*6.242*10**18}
Alpha: {Aa}
Beta: {Ab}

Nun die Intensität an den K-Kanten:
Zink: {I_K_Zink}
Gallium: {I_K_Gallium}
Brom: {I_K_Brom}
Rubidium: {I_K_Rubidium}
Zirkonium: {I_K_Zirkonium}
Strontium: {I_K_Strontium}

Und den Energien:
Zink: {E_K_Zink*6.242*10**18}
Gallium: {E_K_Gallium*6.242*10**18}
Brom: {E_K_Brom*6.242*10**18}
Rubidium: {E_K_Rubidium*6.242*10**18}
Zirkonium: {E_K_Zirkonium*6.242*10**18}
Strontium: {E_K_Strontium*6.242*10**18}

Die Abschirmkonstanten:
Zink: {sigma_K_Zink}
Gallium: {sigma_K_Gallium}
Brom: {sigma_K_Brom}
Rubidium: {sigma_K_Rubidium}
Zirkonium: {sigma_K_Zirkonium}
Strontium: {sigma_K_Strontium}

Für den Rydberg graph und die Ausgleichskurve sind:
m : {unparams[0]}
b : {unparams[1]}
Rydberg-Frequenz: {Ryd_exp}
Rydberg Energie: {Ryd_exp*h}
Rydberg konstante: {Ryd_exp/c}
""")
