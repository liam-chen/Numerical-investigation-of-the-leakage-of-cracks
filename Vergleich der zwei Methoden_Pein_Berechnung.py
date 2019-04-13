# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 15:27:29 2019

@author: LiamChen
"""

from iapws import IAPWS97
import math 
import numpy as np
import matplotlib.pyplot as plt
import csv

def calc_D_h(_A,_U):
    _D_h=4*_A/_U
    print('D_h:',_D_h)
    return _D_h  

def calc_zeta( s, D_h, Rz, zeta_Ein=0.5, zeta_Aus=0): 
   # zeta_Ein : Einlaufverlust 
   # zeta_Aus : Auslaufverlust 
   # _lambda : Widerstandsbeiwert
   #_lambda=1/(max(2,(2**0.5)*math.log10(1.5*D_h*1000/2/Rz)))**2
   _lambda=min(2,1/(2*(math.log10(1.5*D_h*1000/2/Rz))**2))
   zeta = zeta_Ein + _lambda * s / D_h + zeta_Aus
   print('lambda:',_lambda,'zeta',zeta)
   return zeta

def calc_G_HEM_Zweiphasig(p_0, T_0, p_2, T_2, zeta_):  #zweiphasig
    schritt=300    #wie viel Teile geteilt sind in der Gleichung B3.2-9
    p = np.linspace(p_0,p_2,schritt)
    delta_p = (p_2-p_0)/schritt   
    state_0 = IAPWS97(T = T_0, P = p_0)   
    s_0 = state_0.s# unterschiedlicher Zustand mit gleicher Entropie
    sigma_p=0
    sigma_v=0
    G_alt=0
    G_neu=0
    kritisch=False
    for pressureVAL in p: 
        state_1 = IAPWS97(P = pressureVAL, s=s_0)    # unterschiedlicher Zustand mit gleicher Entropie    
        state_2 = IAPWS97(P = (pressureVAL+delta_p), s=s_0)  # unterschiedlicher Zustand mit gleicher Entropie
        v_1= state_1.v
        v_2= state_2.v          
        delta_v=v_2-v_1
        sigma_v=sigma_v+delta_v/v_2
        sigma_p=sigma_p+delta_p/v_2
#        print((-2*sigma_p)/(zeta+2*sigma_v))
        G_neu=(((-2*sigma_p)/(zeta_+2*sigma_v))**0.5)*1000  #kg/m^2/s
        if G_neu>G_alt:
            G_alt=G_neu
        else:
            kritisch=True
            break
    return [round(G_alt,2),round(pressureVAL,4),kritisch]

def calc_P_ein_mit_gleich_T(p_0, T_0, p_2, T_2, zeta_,A_): #beruecksichtigung des Eintrittsducks mit klassischen Formulierung
    state_s_T_0 = IAPWS97(T = T_0, x = 0) 
    state_0 = IAPWS97(T = T_0, P = p_0)
    p_s_0=state_s_T_0.P 
    rho_0=state_0.rho  #dichte kg/m^3
    P_ein_List=np.linspace(max(p_s_0,p_0-1,p_2),p_0,20)
    delta_G_List=[]
    for P_ein in P_ein_List:
        G_1=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)[0]  
        G_2=math.sqrt(2*(p_0-P_ein)*1000000*rho_0)
        print(P_ein,G_1,G_2)
        delta_G=abs(G_1-G_2)
        delta_G_List.append(delta_G)
    delta_G_min=min(delta_G_List)
    G_index=delta_G_List.index(delta_G_min)
    P_ein=P_ein_List[G_index] 
    Zustand=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)
    G=Zustand[0]
    P_aus=Zustand[1]
    kritisch=Zustand[2]
    V=G*A_*1000  #Massenstrom g/s
    print(p_0,round(G,2),round(V,2),P_ein,round(P_aus,4),kritisch)
    return [round(G,2),round(V,2),P_ein,round(P_aus,4),kritisch]

def calc_P_ein_mit_gleich_s(p_0, T_0, p_2, T_2, zeta_,A_): #beruecksichtigung des Eintrittsducks mit der Erweiterungen der klassischen Formulierung
    state_s_T_0 = IAPWS97(T = T_0, x = 0) 
    state_0 = IAPWS97(T = T_0, P = p_0)
    h_0 = state_0.h
    s_0 = state_0.s#
    p_s_0=state_s_T_0.P 
#    rho_0=state_0.rho  #dichte kg/m^3
    P_ein_List=np.linspace(max(p_s_0,p_0-1,p_2),p_0,20)
    delta_G_List=[]
    for P_ein in P_ein_List:
        G_1=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)[0]
        state_1= IAPWS97(P=P_ein,s=s_0)
        h_1=state_1.h
        v_1=state_1.v
        G_2=math.sqrt(2*(h_0-h_1)*1000/v_1/v_1)
        print(P_ein,G_1,G_2)
        delta_G=abs(G_1-G_2)
        delta_G_List.append(delta_G)
    delta_G_min=min(delta_G_List)
    G_index=delta_G_List.index(delta_G_min)
    P_ein=P_ein_List[G_index] 
    Zustand=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)
    G=Zustand[0]
    P_aus=Zustand[1]
    kritisch=Zustand[2]
    V=G*A_*1000  #Massenstrom g/s
    print(p_0,round(G,2),round(V,2),P_ein,round(P_aus,4),kritisch)
    return [round(G,2),round(V,2),P_ein,round(P_aus,4),kritisch]

T_0 = 240 + 273.15  #Stagnationstemperatur  K
p_0 = 4           #Stagnationsdruck   MPa
T_2 = 20 + 273.1    #Umgebungstemperatur  K
p_2 = 0.1           #Umgebungsdruck  MPa

#Geometrie
s=30   #Wanddicke  mm
COD=0.1635  #Rissweite  mm
FCL=27.95  #Risslänge  mm
A=COD*FCL  #Rissfläche mm^2
print('A:%s'%A)
#U=4*math.sqrt((COD/2)**2+(FCL/2)**2)  #Rissumfangslänge  mm  
U=2*(FCL+COD)
D_h=calc_D_h(A,U)  #hydraulischer Durchmesser  mm

#Rz=10  #Rauheit  um
Rz=21.594
zeta=calc_zeta(s,D_h,Rz)   #Stroemungswiderstand  

P_1=calc_P_ein_mit_gleich_T(p_0,T_0,p_2,T_2,zeta,A)
P_2=calc_P_ein_mit_gleich_s(p_0,T_0,p_2,T_2,zeta,A)
print('/n/n')
print(P_1)
print(P_2)