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

def f(x):
    return x**2

x = np.linspace(-10,10,1000)
plt.plot(x,f(x),'r-',label = 'Messdaten')
plt.ylabel(f"f(x)= x*x")
plt.xlabel(f"x")
plt.legend()
plt.savefig("plots/plot.pdf",bbox_inches='tight')
plt.close()


print(f"""


""")
