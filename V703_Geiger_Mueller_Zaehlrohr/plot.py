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

def f(x,m,b):
    return m*x+b

e = scipy.constants.e
U_Kenn,N = np.genfromtxt("data/Kennlinie.dat", unpack = True)
N_Zaehl = np.array([N[3],N[8],N[13],N[18],N[23],N[28],N[33],N[38]])

print("Kennlinie")
N = unp.uarray(N,np.sqrt(N))

Plateau = noms(N[5:-8])
params,cov = curve_fit(f,U_Kenn[5:-8], Plateau)
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)
plt.errorbar(U_Kenn, noms(N), yerr=stds(N), xerr=None,fmt = 'rx', label = "Messwerte")
plt.plot(U_Kenn[5:-8],f(U_Kenn[5:-8],*params),label = "Linearer Fit")
plt.ylabel(f"Impulse pro 60s")
plt.xlabel(f"Spannung in Volt")
plt.legend()
plt.savefig("plots/Kennlinie.pdf",bbox_inches='tight')
plt.close()

print("Totzeit")
N1 = 96041/120
N1 = unp.uarray(N1,np.sqrt(N1))
N12 = 158479/120
N12 = unp.uarray(N12,np.sqrt(N12))
N2 = 76518/120
N2 = unp.uarray(N2,np.sqrt(N2))

Tot = (N1+N2-N12)/(2*N1*N2)

print("Ladungsmenge")
U_Zaehl,I = np.genfromtxt("data/Zaehlrohrstrom.dat", unpack = True)
I = unp.uarray(I,0.05)
I = I*10**(-6)
N_Zaehl = unp.uarray(N_Zaehl,np.sqrt(N_Zaehl))
N_Zaehl = N_Zaehl/60
Z = I/(e*N_Zaehl)

# plt.errorbar(noms(I), noms(Z), yerr=stds(Z), xerr=None,fmt = 'rx', label = "Messwerte")
plt.plot(noms(I), noms(Z),'rx', label = "Messwerte")
plt.ylabel(f"Z")
plt.xlabel(f"I in Ampere")
plt.legend()
plt.savefig("plots/Ladungsmenge.pdf",bbox_inches='tight')
plt.close()

print(f"""
Parameter des Ausgleichs:
m: {unparams[0]}
b: {unparams[1]}

Totzeit:
{Tot}

Ladungen pro Teilchen:
{Z[0]}
{Z[1]}
{Z[2]}
{Z[3]}
{Z[4]}
{Z[5]}
{Z[6]}
{Z[7]}
""")

