# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 12:18:55 2016

@author: Nicolás Torasso - nicolas.torasso@gmail.com

Este es un código para visualizar la solución completa del oscilador
armónico amortiguado y forzado. Se compara con la solución analítica.
Tomo masa = 1
"""
# Importo los módulos a usar de las bibliotecas
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

plt.figure(2)# Abro la ventana de gráfico que voy a usar
plt.clf()#Limpia todas las ventanas de gráfico abiertas

#%% Parámetros relevantes
# Defino los parámetros del FORZANTE: frecuencia, amplitud y fase
w = 1
F = 1
delta = 0
def forzante(t):
    return F*np.cos(w*t+delta)
    
# Defino los parámetros del OSCILADOR: amortiguación y frecuencia propia
gamma = 0.4
w0 = 4

# Defino el vector con las condiciones iniciales (x0 y v0)
X0 = np.array([2,0])    
#%% Solución numérica
# Defino el vector de tiempos para los cuales obtengo la sc numérica
def ecdif(X,t):
    return np.array([X[1] , -gamma*X[1]-w0**2*X[0]+forzante(t)])

t = np.linspace(0,60,10000)
solucion = odeint(ecdif,X0,t)

#%% Solución analítica
# Defino la sc analítica completa
def sc_analitica(t):
    wg = np.sqrt( w0**2-(gamma/2)**2 )
    A = F/np.sqrt( (gamma*w)**2 + (w0**2-w**2)**2 )
    theta = delta-np.arctan(gamma*w/(w0**2-w**2))
    
    c1 = X0[0] - A*np.cos(theta)
    c2 = (X0[1]+gamma*X0[0]/2-gamma*np.cos(theta)/2+A*w*np.sin(theta))/wg

    xh = np.exp(-gamma*t/2)*(c1*np.cos(wg*t) + c2*np.sin(wg*t))    
    xp = A*np.cos(w*t+theta)
    return xh + xp
#%% Gráficos

#plt.plot(t,solucion[:,0],label='Sc. numérica')
plt.plot(t,forzante(t),'k--',label='Forzante')
plt.plot(t,sc_analitica(t),'r-',label='Sc. analítica')

# Agrego nombres al gráficoy a los ejes
plt.xlabel('Tiempo',fontsize=16)
plt.ylabel('Amplitud',fontsize=16)
plt.title('Oscilador forzado',fontsize=16)
plt.legend(loc=0, fontsize=16)
plt.tick_params(labelsize=14)
plt.grid()

plt.show()
