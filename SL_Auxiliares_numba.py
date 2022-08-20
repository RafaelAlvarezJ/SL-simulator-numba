# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 15:29:28 2021

@author: rafa_

Archivo en el cual se encuentran todas las funciones auxiliares para la 
simulación numérica de Streamlines. De esta manera, el código principal no es 
tan extenso."""

import numpy as np
#import SL_Datos as d
from numba import jit, njit, vectorize, prange
import pandas as pd

def LecturaVelocidades():
    Velocidades = pd.read_csv('datos200x200.txt')
    
    U_l = np.array(Velocidades["ux_l"])
    U_r = np.array(Velocidades["ux_r"])
    U_t = np.array(Velocidades["uy_t"])
    U_b = np.array(Velocidades["uy_b"])
    
    aux = len(U_l)
    aux = aux**(1/2)
    aux = int(aux) # aux = nodosx o nodosy, ya que es una malla cuadrada
    
    (U_l, U_r, U_t, U_b) = reacomodo(U_l, U_r, U_t, U_b, aux)
    
    aux2 = aux
    
          #(nodosx, )  
    return (aux, aux2, U_l, U_r, U_t, U_b)
    
#@njit
def reacomodo(U_l, U_r, U_t, U_b, aux):
    U_l = U_l.reshape(aux,aux)
    U_r = U_r.reshape(aux,aux)
    U_t = U_t.reshape(aux,aux)
    U_b = U_b.reshape(aux,aux)
    
    return (U_l, U_r, U_t, U_b)
    

#@jit(nopython=True)
#@njit(parallel=True)
def InitSL(dx,dy,n):
    # Función que nos da los puntos iniciales de las SL, dando como input 
    # el número de SL deseado (debe ser par) y las dimensiones de la celda
    
    # pensada para un modelo five-spot, donde de la celda del pozo inyector, 
    # solo dos de sus caras llevarán streamlines
    
    # - Se asume una celda cuadrada
    
    # n : número de streamlines (tiene que ser un número par)
    
    Xg=np.zeros(n) 
    Yg=np.zeros(n)
    
    aux1 = (n/2) # para los puntos de partida de las SL
    aux2 = (dx/aux1) # para dividir una cara en (aux1) partes iguales
    aux3 = aux2/2
    
    for i in range (n):
        if i < aux1: # para la primera mitad
            Xg[i] = dx
            Yg[i] = aux2*(i)+aux3
        else:        # para la segunda mitad
            Xg[i] = aux2*(i)+aux3-dx
            Yg[i] = dy
    
    return (Xg,Yg)


@njit(cache=True,nogil=True)#,
#@vectorize
def pseudotiempo(c,Dtau):
    "  Función que se ocupa para el calculo de la posición de la partícula "
    if c==0:
        y=Dtau
    else:
        y=(np.exp(c*Dtau)-1)/c
    return (y)

@njit(nogil=True,cache=True)#(parallel=True)
#@vectorize
def velocidad(vp1,cp,p1,p):
    # vp1=velocidad a la iquierda o al fondo de la celda
    # cp=coeficiente
    # p1=coordenada a la izquierda o al fondo de la celda
    # p=punto de partida, (x/y)
    v=vp1+cp*(p-p1)
    return (v)

@njit#(nogil=True)#,cache=True)# y nogil alentan
def cara(csp,i,j):
    """ Función que devuelve la cara entrada actual (cea) de la partícula, también se
    encarga de modificar los índices de celda (i,j) para el siguiente cálculo """ 
    
    # 1 = cara oeste
    # 2 = cara este
    # 3 = cara sur
    # 4 = cara norte
    
    #la cara de salida pasada (csp) es la cara de entrada actual (cea)
    if   csp==1: #Sale al oeste
        cea=2
        i=i-1
        
    elif csp==2: #Sale al este
        cea=1
        i=i+1
        
    elif csp==3: #Sale al sur
        cea=4
        j=j-1
        
    elif csp==4: #Sale al norte
        cea=3
        j=j+1
        
    return (cea,i,j)

#@jit(nopython=True)#, parallel=True) #parallel lo alenta
@njit(cache=True,parallel=True)#nogil=True)#)##)
def algoritmoPollock(x0,y0,ux1,ux2,uy1,uy2,cea,i,j,DX,DY,phi):
    #Cálculo de los coeficientes
    cx=(ux2-ux1)/DX
    cy=(uy2-uy1)/DY
    
    #print(f"i={i}, j={j}")
    #print(f"cx={cx}, cy={cy}, suma={cx+cy}")
    
    
    # cea (cara de entrada actual)
    # 1 = cara oeste
    # 2 = cara este
    # 3 = cara sur
    # 4 = cara norte
    """
    INICIO MODIFICACIÓN
    """
    "La velocidad en x,y es constante."
    if cx==0 and cy==0: 
        if cea==1:   # Cara oeste
            ux0=ux1
            uy0=uy1#velocidad(uy1, cy, j*DY, y0)
            
        elif cea==2: # Cara este
            ux0=ux2
            uy0=uy1#velocidad(uy1, cy, j*DY, y0)
        
        elif cea==3: # Cara sur
            uy0=uy1
            ux0=ux1#velocidad(ux1, cx, i*DX, x0)
            
        elif cea==4: # Cara norte
            uy0=uy2
            ux0=ux1#velocidad(ux1, cx, i*DX, x0)
            
        Dtau_x1=(DX*i    -x0)/ux0    
        Dtau_x2=(DX*(i+1)-x0)/ux0
        Dtau_y1=(DY*j    -y0)/uy0
        Dtau_y2=(DY*(j+1)-y0)/uy0
    
    else: #Si cx,cy son diferentes de 0
            # Cara de entrada actual
        if cea==1:   # Cara oeste
            ux0=ux1
            uy0=velocidad(uy1, cy, j*DY, y0)
            
        elif cea==2: # Cara este
            ux0=ux2
            uy0=velocidad(uy1, cy, j*DY, y0)
        
        elif cea==3: # Cara sur
            uy0=uy1
            ux0=velocidad(ux1, cx, i*DX, x0)
            
        elif cea==4: # Cara norte
            uy0=uy2
            ux0=velocidad(ux1, cx, i*DX, x0)
        
        "Cálculo de tiempo de vuelo a cada cara"
        # Para velocidades en x
        if ux0 == 0: # Para evitar la división entre cero
            Dtau_x1 = np.NaN
            Dtau_x2 = np.NaN
        else:
            if ux1 == 0: # Para evitar el log(cero)
                Dtau_x1 = np.NaN
            else:
                Dtau_x1=(1/cx)*np.log(ux1/ux0)
            
            if ux2 == 0: # Para evitar el log(cero)
                Dtau_x2 = np.NaN
            else:
                Dtau_x2=(1/cx)*np.log(ux2/ux0)        
        
        # Para velocidades en y
        if uy0 == 0: # Para evitar la división entre cero
            Dtau_y1 = np.NaN
            Dtau_y2 = np.NaN
        else:
            if uy1 == 0: # Para evitar el log(cero)
                Dtau_y1 = np.NaN
            else:
                Dtau_y1=(1/cy)*np.log(uy1/uy0)
            
            if uy2 == 0: # Para evitar el log(cero)
                Dtau_y2 = np.NaN 
            else:
                Dtau_y2=(1/cy)*np.log(uy2/uy0)
        
    "Aquí se sale de los condicionales de cx,cy"    
    #print(f"ux1={ux1}, ux2={ux2}, uy1={uy1}, uy2={uy2} ")

    Dtau_array=np.array([Dtau_x1,Dtau_x2,Dtau_y1,Dtau_y2])
    #print(f"{Dtau_array}")

    """
    FIN MODIFICACIÓN
    """

    for i in prange (len(Dtau_array)):
        if Dtau_array[i]<=1:
            Dtau_array[i]=np.NaN
            
    
    Dtau=np.nanmin(Dtau_array)#*phi
    #print (Dtau_array)
    
    #Para devolver la cara de salida 
    if Dtau==Dtau_array[0]:
        csp=1
    elif Dtau==Dtau_array[1]:
        csp=2
    elif Dtau==Dtau_array[2]:
        csp=3
    elif Dtau==Dtau_array[3]:
        csp=4
        
    x_new=x0+ux0*pseudotiempo(cx,Dtau)
    y_new=y0+uy0*pseudotiempo(cy,Dtau)
    
    #print(f"x_new={x_new}, y_new={y_new}")
    #print(f"cea={cea}, csp={csp} \n")
    
    return (x_new,y_new,Dtau*phi, csp)
#    return(Dtau)
#    return(Dtau_array)

@njit(cache=True)#(nogil=True)#cache=True)#,)
def kr(Sw, Srw, Sro):
    omega = 2
    
    "Función que calcula las permeabilidades relativas cuadráticas"
    Scf=(Sw-Srw)/(1-Srw-Sro) # saturación efectiva. Ver Diaz-Viera
    
    krw=Scf**omega
    kro=(1-Scf)**omega
    
    return(krw, kro)

@njit(cache=True)#(nogil=True)#cache=True)#, # vectorize y njit dan resultados similares
#@vectorize
def fractionalflow(kro,krw,Muo,Muw):
    "Función que se encarga de calcular el flujo fraccional de agua"
    
    if krw ==0 :
        fw = 0
    
    else:
        fw = 1/(1+(kro*Muw)/(krw*Muo))
    
    return(fw)
    

@njit(cache=True)#(nogil=True)#cache=True)#,nogil=True)#parallel=True)
def transporte(Sw, fw, Srw, Sro, Muo, Muw, DTAU, dt, TiempoTotal):#, fw_avg):
    "Función que calcula la saturación explícitamente"
    # Sw, DTAU, fw son arreglos. Las demás variables son escalares.
    
    # tiempodevuelo=sum(DTAU)
    
    t=0 
    "Ciclo temporal"
    while t <= TiempoTotal:
        #dfw_aux = np.zeros(len(DTAU)) # Tiene un elemento menos
        
        ### DEJEMOS IGUAL
        for i in range (len(DTAU)): # A lo largo de la streamline. Este ciclo podría solo calcular fw sin Sw.

            if i==0: # Condición de inyección
                Sw[i] = (1-Sro)
                
                krw, kro = kr(Sw[i], Srw, Sro) # kro, krw se actualizan en cada ciclo
                fw[i] = fractionalflow(kro, krw, Muo, Muw)
                
            else: # Cálculo IMPES
                if DTAU[i]==0:
                    break
            
                krw, kro = kr(Sw[i], Srw, Sro) # kro, krw se actualizan en cada ciclo
                fw[i] = fractionalflow(kro, krw, Muo, Muw)
            
                aux = (dt/DTAU[i])*(fw[i]-fw[i-1])
                Sw[i] = Sw[i]-aux
                
        t = t+dt

@njit(nogil=True)#cache=True)#,nogil=True)#,parallel=True)
def principal(Xg, Yg, X_coord, Y_coord, Tau_SL, Sw_SL, limit, U_xl, U_xr, U_yb, U_yt, 
              NSL, xcell_inj, ycell_inj, DX, DY, phi,
              Sw_ini, Srw, Sor, Muo, Muw, 
              dt, TiempoTotal):
    for k in range (NSL): # Para todas las  SL
        
        s = 0 #Contador para los loops del ciclo para cada SL
        i = xcell_inj; j = ycell_inj # es la celda de partida de cada SL (0,0) 

        X_coord[k][0] = Xg[k] # Se inicia la trayectoria de la streamline
        Y_coord[k][0] = Yg[k] # Se inicia la trayectoria de la streamline

        #set_num_threads(4)  #Selecionar numero de threads 
        for l in range (1,limit-2): #while (i*j != (nx_cell-1)*(nx_cell-2) ):

            if s==0: # Condicional para asignar cara de entrada, al inicio del algoritmo para una SL
                if   X_coord[k][0] == DX:
                    csp=2
                elif Y_coord[k][0] == DY:
                    csp=4
        
            cea,i,j =cara(csp,i,j) #aquí se obtiene la cara de entrada actual (cea) de la partícula
            aux1,aux2,aux3,csp = algoritmoPollock(X_coord[k][l-1], Y_coord[k][l-1], 
                                                U_xl[i][j], U_xr[i][j], 
                                                U_yb[i][j], U_yt[i][j], 
                                                cea, i,j, DX, DY, phi)
        
            X_coord[k][l] = aux1 #X = np.append(X,aux1) # x_new
            Y_coord[k][l] = aux2 #Y = np.append(Y,aux2) # y_new
            Tau_SL [k][l] = aux3 #DTAU = np.append(DTAU,aux3) # Dtau*phi
        
            s=s+1
            if s>500: #Máximo de loops por SL 
                break
        
        "Aquí se hacen los cálculos de transporte por streamline"
        DTAU = Tau_SL[k][:] # La lista se hace arreglo para los cálculos de transporte
        fw = np.zeros(len(DTAU))
        Sw = np.ones (len(DTAU))*Sw_ini
        transporte(Sw, fw, Srw, Sor, Muo, Muw, DTAU, dt, TiempoTotal)#, fw_avg)
        Sw_SL[k][:]=Sw #Sw_SL.append(Sw)

