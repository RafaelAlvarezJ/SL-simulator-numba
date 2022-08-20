# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 15:43:04 2021

@author: rafa_

Archivo en el que se ingresan los datos para la simulación numérica de Streamlines.
"""
import Units as u

NSL=40 # Número de Streamlines

"Dimensiones de la malla"
nx_cell=200 # celdas en x
ny_cell=200 # celdas en y

"Datos de simulación"
dt = 0.001*u.day
TiempoTotal = 100*u.day #100*u.day

"Propiedades del medio"
LX = 1000*u.ft #ft
DX = LX/nx_cell
#x=np.linspace(0,250,6)
LY = 1000*u.ft #ft
DY = LY/ny_cell
#y=np.linspace(0,250,6)
DZ = 100*u.ft #ft

phi = 0.2
Sw_ini = 0.22

"Propiedades de los fluidos"
Muo = 1.14*u.Cp
Muw = 0.096*u.Cp

Swc = 0.2
Srw = 0.22
Sor = 0.2

"Pozos"
Q_inj  = 50 #bl/d
Q_prod = 100 #bl/d

# Celdas de los pozos (índices)
xcell_inj  = 0; ycell_inj  = 0 # Inyección
xcell_prod = nx_cell-1; ycell_prod = ny_cell-1 # Producción

"Matrices de Velocidades, obtenidas de un simulador de diferencias finitas"
#Left
# U_xl=np.array([[0,0,0,0,0],
#                 [0.28072,0.11485,0.0638 ,0.05104,0.05104],
#                 [0.16589,0.11485,0.08932,0.08932,0.10209],
#                 [0.10209,0.08932,0.08932,0.11485,0.16588],
#                 [0.05104,0.05104,0.06381,0.11484,0.28073]])


#Right
# U_xr=np.array([[0.28072,0.11485,0.0638 ,0.05104,0.05104],
#                 [0.16589,0.11485,0.08932,0.08932,0.10209],
#                 [0.10209,0.08932,0.08932,0.11485,0.16588],
#                 [0.05104,0.05104,0.06381,0.11484,0.28073],
#                 [0,0,0,0,0]])

#Bottom
# U_yb=np.array([[0,0.28072,0.16589,0.10209,0.05104],
#                 [0,0.11485,0.11485,0.08932,0.05104],
#                 [0,0.0638 ,0.08932,0.08932,0.06381],
#                 [0,0.05104,0.08932,0.11485,0.11484],
#                 [0,0.05104,0.10209,0.16588,0.28073]])



#Top
# U_yt=np.array([[0.28072,0.16589,0.10209,0.05104,0],
#                 [0.11485,0.11485,0.08932,0.05104,0],
#                 [0.0638 ,0.08932,0.08932,0.06381,0],
#                 [0.05104,0.08932,0.11485,0.11484,0],
#                 [0.05104,0.10209,0.16588,0.28073,0]])



"Puntos iniciales de las SL"


# Xg=np.array([50 , 50, 50  , 41.7, 25 , 8.3])
# Yg=np.array([8.3, 25, 41.7, 50  , 50 , 50 ])

