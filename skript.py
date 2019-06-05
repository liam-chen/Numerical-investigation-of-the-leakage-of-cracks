# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 09:31:26 2018

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

def calc_G_Benourlli(p_0, T_0, p_2, T_2, zeta_): #einphasig

    state_0 = IAPWS97(T = T_0, P = p_0)

    rho_0=state_0.rho
    G = math.sqrt(2*(p_0-p_2)*1000000*rho_0/(1+zeta_)) 
   
    return round(G,2)

def calc_G_modB_Einphasig(p_0,T_0,zeta_): #einphasig
   state_s_T_0 = IAPWS97(T = T_0, x = 0) 
   p_s_T_0 = state_s_T_0.P 
   rho_s_T_0 = state_s_T_0.rho   #Density, [kg/m³]
   G = math.pow(2 * (p_0 - p_s_T_0)*1000000 * rho_s_T_0/(1+zeta_), 0.5)
   
   return round(G,2)

def calc_G_HEM_Zweiphasig_1(p_0, T_0, p_2, T_2, zeta_):  #kritisch Massenstrom darstellen
    schritt=1000    #wie viel Teile geteilt sind in der Gleichung B3.2-9
    p = np.linspace(p_0,p_2,schritt)
    delta_p = (p_2-p_0)/schritt   
    state_0 = IAPWS97(T = T_0, P = p_0)   
    s_0 = state_0.s # unterschiedlicher Zustand mit gleicher Entropie
    v_0 = state_0.v
    sigma_p=0
    sigma_v=0
    G=0
    G_List=[]
    for pressureVAL in p: 
#        state_1 = IAPWS97(P = pressureVAL, s=s_0)    # unterschiedlicher Zustand mit gleicher Entropie    
        state_2 = IAPWS97(P = (pressureVAL+delta_p), s=s_0)  # unterschiedlicher Zustand mit gleicher Entropie
#        v_1= state_1.v
        v_2= state_2.v          
#        delta_v=v_2-v_1
        sigma_v=math.log(v_2/v_0)
        sigma_p=sigma_p+delta_p/v_2
        G=(((-2*sigma_p)/(zeta_+2*sigma_v))**0.5)*1000  #kg/m^2/s
        G_List.append(G)
    plt.figure(1)   #Diagramm zeichnen
    plt.xlabel("Austrittsdruck in MPa")
    plt.ylabel("G")
    plt.plot(p, G_List)
    return [round(G,2),round(pressureVAL,4)]

def calc_G_HEM_Zweiphasig(p_0, T_0, p_2, T_2, zeta_):  #zweiphasig    genauig Ergebnisse
    schritt=1000    #wie viel Teile geteilt sind in der Gleichung B3.2-9
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
    return [round(G_alt,2),round(pressureVAL,4),kritisch]


def calc_P_ein_und_G_HEM(p_0, T_0, p_2, T_2, zeta_,A_): #beruecksichtigung des Eintrittsducks  mit Bisektionsverfahren
    state_s_T_0 = IAPWS97(T = T_0, x = 0) 
    state_0 = IAPWS97(T = T_0, P = p_0)
    p_s_0=state_s_T_0.P 
    rho_0=state_0.rho  #dichte kg/m^3
    a=max(p_s_0,p_0-1,p_2)
    b=p_0
    delta_G=calc_G_HEM_Zweiphasig(a,T_0,p_2,T_2,zeta_)[0]- math.sqrt(2*(p_0-a)*1000000*rho_0) 
    P_ein_List=[]
    G1_List=[]
    G2_list=[]
    while abs(delta_G)>10:
        c=(a+b)/2
        delta_G=calc_G_HEM_Zweiphasig(c,T_0,p_2,T_2,zeta_)[0]- math.sqrt(2*(p_0-c)*1000000*rho_0) 
        if  delta_G>0:
            b=c
        else:
            a=c
        print(c,math.sqrt(2*(p_0-c)*1000000*rho_0) )
        P_ein_List.append(c)
        G1_List.append(calc_G_HEM_Zweiphasig(c,T_0,p_2,T_2,zeta_)[0])
        G2_list.append(math.sqrt(2*(p_0-c)*1000000*rho_0))
    P_ein=c
    Zustand=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)
    G=Zustand[0]
    P_aus=Zustand[1]
    kritisch=Zustand[2]
    V=G*A_/1000000  #Massenstrom kg/m^2/s
    state_2=IAPWS97(T = T_2, P = p_2)
    rho_2=state_2.rho
#    print(rho_2)
    v=G/rho_2   # m/s
    #print(p_0,round(G,2),round(V,6),P_ein,round(P_aus,4),kritisch)
    #print('v(m/s)=',v)    #durchschnittliche Geschwindigkeit
#    print(P_ein_List, G1_List, G2_list)
#    plt.figure(2)   #Diagramm zeichnen
#    plt.xlabel("Eintrittsdruck in MPa")
#    plt.ylabel("G")
#    plt.plot(P_ein_List, G1_List)  
#    plt.plot(P_ein_List, G2_list)   
#    
    
    return [round(G,2),round(V,6),P_ein,round(P_aus,4),kritisch] #G:Leckrate  V:Massenstrom


#def calc_P_ein_und_G_HEM_H(p_0, T_0, p_2, T_2, zeta_,A_): #beruecksichtigung des Eintrittsducks  mit Bisektionsverfahren
#    state_s_T_0 = IAPWS97(T = T_0, x = 0) 
#    state_0 = IAPWS97(T = T_0, P = p_0)
#    h_0 = state_0.h
#    s_0 = state_0.s
#    p_s_0=state_s_T_0.P 
#    rho_0=state_0.rho  #dichte kg/m^3
#    a=max(p_s_0,p_0-1,p_2)
#    b=p_0
#    state_a = IAPWS97( s=s_0  ,P=a )
#    h_a = state_a.h
#    v_a = state_a.v
#    delta_G=calc_G_HEM_Zweiphasig(a,T_0,p_2,T_2,zeta_)[0]- math.sqrt(2*(h_0-h_a)*1000/v_a/v_a)
#    P_ein_List=[]
#    G1_List=[]
#    G2_list=[]
#    while abs(delta_G)>10:
#        c=(a+b)/2
#        state_c = IAPWS97( s=s_0  ,P=c )
#        h_c = state_c.h
#        v_c = state_c.v
#        delta_G=calc_G_HEM_Zweiphasig(c,T_0,p_2,T_2,zeta_)[0]- math.sqrt(2*(h_0-h_c)*1000/v_c/v_c)
#        if  delta_G>0:
#            b=c
#        else:
#            a=c
#        print(c,math.sqrt(2*(h_0-h_c)*1000/v_c/v_c))
#        P_ein_List.append(c)
#        G1_List.append(calc_G_HEM_Zweiphasig(c,T_0,p_2,T_2,zeta_)[0])
#        G2_list.append(math.sqrt(2*(p_0-c)*1000000*rho_0))
#    P_ein=c
#    Zustand=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)
#    G=Zustand[0]
#    P_aus=Zustand[1]
#    kritisch=Zustand[2]
#    V=G*A_/1000000  #Massenstrom kg/m^2/s
#    state_2=IAPWS97(T = T_2, P = p_2)
#    rho_2=state_2.rho
##    print(rho_2)
#    v=G/rho_2   # m/s
#    #print(p_0,round(G,2),round(V,6),P_ein,round(P_aus,4),kritisch)
#    #print('v(m/s)=',v)    #durchschnittliche Geschwindigkeit
##    print(P_ein_List, G1_List, G2_list)
##    plt.figure(2)   #Diagramm zeichnen
##    plt.xlabel("Eintrittsdruck in MPa")
##    plt.ylabel("G")
##    plt.plot(P_ein_List, G1_List)  
##    plt.plot(P_ein_List, G2_list)   
##    
#    
#    return [round(G,2),round(V,6),P_ein,round(P_aus,4),kritisch] #G:Leckrate  V:Massenstrom


T_0 = 101 + 273.15  #Stagnationstemperatur  K
p_0 = 5          #Stagnationsdruck   MPa
T_2 = 20 + 273.1    #Umgebungstemperatur  K
p_2 = 0.1           #Umgebungsdruck  MPa

#Geometrie
s=30   #Wanddicke  mm
COD=0.1635  #Rissweite  mm
FCL=27.95 #Risslänge  mm
A=COD*FCL  #Rissfläche mm^2
print('A:%s'%A)
#U=4*math.sqrt((COD/2)**2+(FCL/2)**2)  #Rissumfangslänge  mm  
U=2*(FCL+COD)
D_h=calc_D_h(A,U)  #hydraulischer Durchmesser  mm

#Rz=10  #Rauheit  um
Rz=21.594
zeta=calc_zeta(s,D_h,Rz)   #Stroemungswiderstand  

#calc_G_HEM_Zweiphasig_1(p_0, T_0, p_2, T_2, zeta)
print('T_0=',T_0-273.15,'\tp_0=',p_0)
print('HEM',calc_P_ein_und_G_HEM(p_0, T_0, p_2, T_2, zeta,A))
###print('HEM_H',calc_P_ein_und_G_HEM_H(p_0, T_0, p_2, T_2, zeta,A))
print('Beourlli',calc_G_Benourlli(p_0, T_0, p_2, T_2, zeta))
print('Modifiziert B',calc_G_modB_Einphasig(p_0, T_0, zeta))