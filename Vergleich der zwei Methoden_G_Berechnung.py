# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 15:03:53 2019

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
    schritt=500    #wie viel Teile geteilt sind in der Gleichung B3.2-9
    p = np.linspace(p_0,p_2,schritt)
    delta_p = (p_2-p_0)/schritt   
    state_0 = IAPWS97(T = T_0, P = p_0)   
    s_0 = state_0.s # unterschiedlicher Zustand mit gleicher Entropie
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
        G_neu=((-2*sigma_p)/(zeta+2*sigma_v))**0.5*1000  #kg/m^2/s
        if G_neu>G_alt:
            G_alt=G_neu
        else:
            kritisch=True
            break
    return [G_alt,pressureVAL,kritisch]


def calc_G_HEM_Zweiphasig_2(p_0, T_0, p_2, T_2, zeta_):  #zweiphasig
    schritt=500    #wie viel Teile geteilt sind in der Gleichung B3.2-9
    p = np.linspace(p_0,p_2,schritt)
    delta_p = (p_2-p_0)/schritt   
    state_0 = IAPWS97(T = T_0, P = p_0)   
    s_0 = state_0.s # unterschiedlicher Zustand mit gleicher Entropie
    v_0 = state_0.v
    sigma_p=0
    sigma_v=0
    G_alt=0
    G_neu=0
    kritisch=False
    for pressureVAL in p: 
#        state_1 = IAPWS97(P = pressureVAL, s=s_0)    # unterschiedlicher Zustand mit gleicher Entropie    
        state_2 = IAPWS97(P = (pressureVAL+delta_p), s=s_0)  # unterschiedlicher Zustand mit gleicher Entropie
#        v_1= state_1.v
        v_2= state_2.v          
#        delta_v=v_2-v_1
        sigma_v=math.log(v_2/v_0)
        sigma_p=sigma_p+delta_p/v_2
        G_neu=((-2*sigma_p)/(zeta+2*sigma_v))**0.5*1000  #kg/m^2/s
        if G_neu>G_alt:
            G_alt=G_neu
        else:
            kritisch=True
            break
    return [G_alt,pressureVAL,kritisch]



T_0 = 300 + 273.15  #Stagnationstemperatur  K
p_0 = 15           #Stagnationsdruck   MPa
T_2 = 20 + 273.1    #Umgebungstemperatur  K
p_2 = 0.1           #Umgebungsdruck  MPa
state_0 = IAPWS97(T = T_0, P = p_0)
rho_0=state_0.rho  #dichte kg/m^3
state_s_T_0 = IAPWS97(T = T_0, x = 0) 
p_s_0=state_s_T_0.P
print(p_s_0)
#Geometrie
s=50   #Wanddicke  mm
COD=1.6  #Rissweite  mm
FCL=650  #Risslänge  mm
A=COD*FCL/2  #Rissfläche mm^2
print('flaeche:',A)
U=4*math.sqrt((COD/2)**2+(FCL/2)**2)  #Rissumfangslänge  mm  
D_h=calc_D_h(A,U)  #hydraulischer Durchmesser  mm

Rz=10  #Rauheit  um

zeta=calc_zeta(s,D_h,Rz)   #Stroemungswiderstand  


G_1=calc_G_HEM_Zweiphasig(p_0,T_0,p_2,T_2,zeta)
G_2=calc_G_HEM_Zweiphasig_2(p_0,T_0,p_2,T_2,zeta)
print(G_1)
print(G_2)