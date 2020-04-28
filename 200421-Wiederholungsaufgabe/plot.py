import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp 

linie, spannung = np.genfromtxt("data.txt",unpack=True)

D = (linie-1)*6
print(f'Abstand D: {D} in mm')
D = D*0.001
print(f'Abstand D: {D} in m')
print(f'Spannung U: {spannung} in V')
def f(m,x,b):
    return m*x+b

params, cov = curve_fit(f,spannung,D)
errors = np.sqrt(np.diag(cov))
params_error = unp.uarray(params,errors)
print(f'Steigung: {params_error[0]}')
print(f'Achsenabschnitt: {params_error[1]}')

x_plot = np.linspace(spannung[0],spannung[8])
plt.plot(x_plot,f(x_plot,*params))
plt.plot(spannung,D,"rx")
plt.xlabel(f"$U$ in V")
plt.ylabel("$D$ in m")
plt.savefig('build/plot.pdf',bbox_inches='tight')