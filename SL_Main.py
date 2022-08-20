# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 21:50:32 2021

@author: rafa_
"""
#Programa para el ejercicio 3.30 del libro Datta-Gupta

import numpy as np
#import pandas as pd
import Units as u
import SL_Auxiliares_numba as A
import SL_Datos as d
#from numba import njit, set_num_threads


from time import perf_counter

import matplotlib.pyplot as plt

#from matplotlib.collections import LineCollection

#start = perf_counter()
(nx_cell, ny_cell, U_xl, U_xr, U_yt, U_yb) = A.LecturaVelocidades()

limit = 400 # maxima loongitud de una línea de corriente
"Aquí guardaremos las trayectorias y saturaciones de todas las streamlines"
X_coord     = np.zeros([d.NSL,limit]) # Aquí guardaremos las coordenadas x            de todas las SL
Y_coord     = np.zeros([d.NSL,limit]) # Aquí guardaremos las coordenadas y            de todas las SL
Sw_SL       = np.zeros([d.NSL,limit]) # Aquí guardaremos la Saturación                de todas las SL
Tau_SL      = np.zeros([d.NSL,limit]) # Aqui guardaremos el tiempo de vuelo acumulado de todas las SL

# Celdai      = np.array([]) # Aquí se guarda las celda i por donde van pasando las SL
# Celdaj      = np.array([]) # Aquí se guarda las celda j por donde van pasando las SL
# auxTau      = np.array([]) # Aquí se guarda (tau*phi) acumulado           de  las SL
# auxDTau     = np.array([]) # Aquí se guarda (dtau*phi)                    de  las SL
# TauCentrado = np.array([]) # QUITALO

start = perf_counter()
"""
CÁLCULO DE PUNTOS INICIALES DE LAS STREAMLINES
"""
Xg,Yg = A.InitSL(d.DX, d.DY, d.NSL)

"""
INICIO DEL ALGORITMO DE POLLOCK
"""
#Tiempo=np.arange(0,d.TiempoTotal,d.dt)
#fw_avg=[] # Lista para la gráfica de fw
# fw_avg = np.array([])
# TAU=[] # Lista donde se van guardando los diversos dtau de la streamline en cuestión
# ejeXgraficoH=[] # Lista creada para numerar las streamlines en el gráfico de tiempo de vuelo

"""
PRINCIPAL
"""
#X_coord[:][0]=Xg
#Y_coord[:][0]=Yg

    
"""
PRINCIPAL
"""
A.principal(Xg, Yg, X_coord, Y_coord, Tau_SL, Sw_SL, limit, U_xl, U_xr, U_yb, U_yt, 
            d.NSL, d.xcell_inj, d.ycell_inj, d.DX, d.DY, d.phi,
            d.Sw_ini, d.Srw, d.Sor, d.Muo, d.Muw, 
            d.dt, d.TiempoTotal)


end = perf_counter()
print(f"Tiempo de ejecución: {end-start} [s]")

"Post proceso sin txt"
# Sw_SL2      = np.asarray(Sw_SL) ; Sw_SL2 = Sw_SL2.flatten()

# DATOS = pd.DataFrame(Celdai)#, columns="Celdai") 
# DATOS["CeldaY"] = pd.DataFrame(Celdaj)
# DATOS["Tau"] = pd.DataFrame(auxTau)
# DATOS["Dtau"] = pd.DataFrame(auxDTau)
# DATOS["TauCentr"] = pd.DataFrame(TauCentrado)
# DATOS["Sw"] = pd.DataFrame(Sw_SL2)

# Columnas = ["CeldaX","CeldaY","Tau","Dtau","TauCentr","Sw"] # Del dataframe
# DATOS.columns = Columnas

# DATOS["Sw*Dtau"]=DATOS["Sw"]*DATOS["Dtau"]

# # Se crean los arrays para los gráficos de contorno
# graftau = np.zeros([nx_cell, ny_cell])
# grafSw  = np.zeros([nx_cell, ny_cell])

# "INICIO MODIFICACIÓN"
# for i in range (nx_cell):
#     for j in range (ny_cell):
        
#         if (i==0) and (j==0): # en el pozo inyector
#             graftau[i][j] = 0
#             grafSw[i][j]  = (1-d.Sor) 
            
#         elif (i == nx_cell-1) and (j == ny_cell-1): # en el pozo productor
            
#             # Para el tiempo de vuelo
#             "Opción 1"
#             array_TAU=np.array(TAU)
#             graftau[i][j] = array_TAU.mean()
            
#             "Opción 2"
#             # graftau[i][j] = DATOS["Tau"].max()
            
#             # Para la saturación
#             grafSw[i][j] = DATOS["Sw"].min()
        
#         else:
            
#             df_aux = DATOS[ ( DATOS["CeldaX"].isin([i]) ) & ( DATOS["CeldaY"].isin([j]) ) ]
            
#             # Para el tiempo de vuelo
#             graftau[i][j] = df_aux["TauCentr"].mean()
            
#             # Para la saturación
#             numerador = df_aux["Sw*Dtau"].sum()
#             denominador = df_aux["Dtau"].sum()
#             grafSw[i][j] = numerador/denominador

"FIN MODIFICACIÓN"

#fw_avg2=A.promedio_fw(fw_avg,Tiempo,d.NSL)

"""
GRÁFICAS
"""

"Gráfica de Streamlines"
#plt.figure("Gráfica de Streamlines", figsize=(10,8))    
# fig, ax = plt.subplots()
# ax.set_xlim(0 , d.LX )
# ax.set_ylim(0 , d.LY ) 


"Gráfica de contorno 1"
# holax=np.linspace(0, d.LX, nx_cell)
# holay=np.linspace(0, d.LY, ny_cell)
# xg,yg = np.meshgrid(holax,holay) # Generación de malla para graficar

# plt.figure("Tiempo de Vuelo", figsize = (10,8))
# plt.contourf(xg,yg,graftau/u.day)
# plt.colorbar()
# plt.title("Tiempo de Vuelo, $[d]$")
# plt.xlabel("$x, [m]$")
# plt.ylabel("$y, [m]$")
# #plt.grid()

# plt.figure("Saturación de agua", figsize = (10,8))
# plt.contourf(xg,yg,grafSw)
# plt.colorbar()
# plt.title("Saturación de agua")
# plt.xlabel("$x, [m]$")
# plt.ylabel("$y, [m]$")
# plt.grid()


"""
Gráficas con matplotlib
"""

plt.figure(f'Saturación de agua a {int(d.TiempoTotal/u.day)} días')#, figsize=(4.8,4.8) )
#plt.figure("SLS_90_3")#, figsize=(4.8,4.8))
for i in range (d.NSL):
    #plt.plot(X_coord[i][:], Y_coord[i][:], c="black",linewidth=0.2)#,marker=".")
    
    sc=plt.scatter(X_coord[i][:], Y_coord[i][:], c=Sw_SL[i][:], s=50, marker="o")
    #plt.tricontourf(X_coord[i][:], Y_coord[i][:], c=Sw_SL[i][:])
    
plt.colorbar(sc)

#plt.title(f'Líneas de Corriente. $S_w$ a {int(d.TiempoTotal/u.day)} [d]')
plt.title(f'Sw a {int(d.TiempoTotal/u.day)} [d]')
#plt.title("Trayectorias")

plt.xlabel("$x, [m]$")
plt.xlim(0,d.LX)
plt.ylabel("$y, [m]$")
plt.ylim(0,d.LY)

#plt.savefig("SLS_90_3.png", dpi=300)
#plt.grid()
