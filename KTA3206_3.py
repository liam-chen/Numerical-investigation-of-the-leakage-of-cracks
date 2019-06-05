# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 13:04:02 2018

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


def calc_G_HEM_Zweiphasig_alt(p_0, T_0, p_2, T_2, zeta_):  #zweiphasig
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

def calc_G_HEM_Zweiphasig(p_0, T_0, p_2, T_2, zeta_):  #zweiphasig    genauiger als frueher 
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

def calc_P_ein_und_G_HEM_alt_nicht_genau(p_0, T_0, p_2, T_2, zeta_,A_): #beruecksichtigung des Eintrittsducks
    state_s_T_0 = IAPWS97(T = T_0, x = 0) 
    state_0 = IAPWS97(T = T_0, P = p_0)
    p_s_0=state_s_T_0.P 
    rho_0=state_0.rho  #dichte kg/m^3
    P_ein_List=np.linspace(max(p_s_0,p_0-0.05,p_2),p_0,21)
    #print(P_ein_List)
    delta_G_List=[]
    for P_ein_variable in P_ein_List:
        G_1=calc_G_HEM_Zweiphasig_2(P_ein_variable,T_0,p_2,T_2,zeta_)[0]  
        G_2=math.sqrt(2*(p_0-P_ein_variable)*1000000*rho_0)
        delta_G=abs(G_1-G_2)
        delta_G_List.append(delta_G)
    delta_G_min=min(delta_G_List)
    G_index=delta_G_List.index(delta_G_min)
    P_ein=P_ein_List[G_index] 
    Zustand=calc_G_HEM_Zweiphasig_2(P_ein,T_0,p_2,T_2,zeta_)
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
    return [round(G,2),round(V,6),P_ein,round(P_aus,4),kritisch] #G:Leckrate  V:Massenstrom

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
        #print(c,math.sqrt(2*(p_0-c)*1000000*rho_0) )
        P_ein_List.append(c)
        G1_List.append(calc_G_HEM_Zweiphasig(c,T_0,p_2,T_2,zeta_)[0])
        G2_list.append(math.sqrt(2*(p_0-c)*1000000*rho_0))
    P_ein=c
    Zustand=calc_G_HEM_Zweiphasig(P_ein,T_0,p_2,T_2,zeta_)
    G=Zustand[0]  #Massenstrom kg/m^2/s
    P_aus=Zustand[1]
    kritisch=Zustand[2]
    V=G*A_/1000   #g/s
    state_2=IAPWS97(T = T_2, P = p_2)
    rho_2=state_2.rho
    print(p_0,T_0,G,V)
    return [round(G,2),P_ein,round(P_aus,4),kritisch,p_s_0] #G:Leckrate  V:Massenstrom


T_0 = 170 + 273.15  #Stagnationstemperatur  K
p_0 = 1          #Stagnationsdruck   MPa
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

#einmal Berechnung
#print(calc_G_modB_Einphasig(13.8,T_0,zeta))
#print(calc_G_HEM_Zweiphasig(p_0,T_0,p_2,T_2,zeta))
#print(calc_G_HEM_Zweiphasig_alt(p_0,T_0,p_2,T_2,zeta))
print(calc_P_ein_und_G_HEM(p_0,T_0,p_2,T_2,zeta,A))
#
#print('\n','modB')
#
#print(calc_G_modB_Einphasig(p_0,T_0,zeta))


##unterschiedlicher Stagnationsdruck mit Beruecksichtigung des Eintrittsdruck
#G_List=[]
#Pstag_List=np.linspace(1.5,6,10  )  #Die Serie des Stagnationsdrucks  Kann einen Fehler melden. falls Pstag ist 
##Pstag_List=[0.2]+list(Pstag_List)
#print(Pstag_List)
#headers = ['Pstag','T','HEM_G','v','Pein','Paus','kritisch','saettigungsdruck']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for P_stag in Pstag_List:
#    Zustand=calc_P_ein_und_G_HEM(P_stag,T_0,p_2,T_2,zeta,A)
#    G_List.append(Zustand[0])
#    row=[P_stag]+[T_0-273.15]+Zustand
#    rows.append(row)
#with open('Dh_HEM_180.csv','w',newline ='') as f:  #schreiben in ein Datei
#    f_csv = csv.writer(f)
#    f_csv.writerow(headers)
#    f_csv.writerows(rows)
#plt.figure(1)   #Diagramm zeichnen
#plt.xlabel("Pstag in MPa")
#plt.ylabel("G")
#plt.plot(Pstag_List, G_List)

##unterschiedlicher Temperatur
#G_List=[]
#Tstag_List=np.linspace(220,277,20)  #Die Serie des Stagnationsdrucks  Kann einen Fehler melden. falls Pstag ist 
##Pstag_List=[1.75]+list(Pstag_List)
#print(Tstag_List)
#headers = ['Pstag','T','Benourlli','modB','HEM','Pein','Paus','kritisch','saettigungsdruck']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for T_stag in Tstag_List:
#    G_b=calc_G_Benourlli(p_0, T_stag+ 273.15, p_2, T_2, zeta)
#    G_m_b=calc_G_modB_Einphasig(p_0,T_stag+ 273.15,zeta)
#    G_HEM=calc_P_ein_und_G_HEM(p_0, T_stag+ 273.15, p_2, T_2, zeta,A)
##    G_List.append(Zustand)
#    row=[p_0]+[T_stag]+[G_b]+[G_m_b]+G_HEM
#    rows.append(row)
#with open('Dh_HEM.csv','w',newline ='') as f:  #schreiben in ein Datei
#    f_csv = csv.writer(f)
#    f_csv.writerow(headers)
#    f_csv.writerows(rows)
##plt.figure(1)   #Diagramm zeichnen
##plt.xlabel("Pstag in MPa")
##plt.ylabel("G")
##plt.plot(Pstag_List, G_List)

##unterschiedlicher Stagnationsdruck mit Bernoulli
#G_List=[]
##Pstag_List=np.linspace(0.5,6,12)  #Die Serie des Stagnationsdrucks  Kann einen Fehler melden. falls Pstag ist 
##Pstag_List=[0.2]+list(Pstag_List)
#print(Pstag_List)
#headers = ['Pstag','G']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for P_stag in Pstag_List:
#    Zustand=calc_G_Benourlli(P_stag,T_0,p_2,T_2,zeta)
#    G_List.append(Zustand)
#    row=[P_stag]+[Zustand]
#    rows.append(row)
#with open('Dh_Benourlli_1.csv','w',newline ='') as f:  #schreiben in ein Datei
#    f_csv = csv.writer(f)
#    f_csv.writerow(headers)
#    f_csv.writerows(rows)
#plt.figure(1)   #Diagramm zeichnen
#plt.xlabel("Pstag in MPa")
#plt.ylabel("G")
#plt.plot(Pstag_List, G_List)
###
##
##unterschiedlicher Stagnationsdruck mit modifizierte Bernoulli
#G_List=[]
##Pstag_List=np.linspace(0.5,6,12)  #Die Serie des Stagnationsdrucks  Kann einen Fehler melden. falls Pstag ist 
##Pstag_List=[0.2]+list(Pstag_List)
#print(Pstag_List)
#headers = ['Pstag','G']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for P_stag in Pstag_List:
#    Zustand=calc_G_modB_Einphasig(P_stag,T_0,zeta)
#    G_List.append(Zustand)
#    row=[P_stag]+[Zustand]
#    rows.append(row)
#with open('Dh_mod_Benourlli_1.csv','w',newline ='') as f:  #schreiben in ein Datei
#    f_csv = csv.writer(f)
#    f_csv.writerow(headers)
#    f_csv.writerows(rows)
#plt.figure(1)   #Diagramm zeichnen
#plt.xlabel("Pstag in MPa")
#plt.ylabel("G")
#plt.plot(Pstag_List, G_List)




##unterschiedlicher Stagnationsdruck ohne berücksichtigung des Eitrittsdrucks
#G_List=[]
#Pstag_List=np.linspace(9,15,21)  #Die Serie des Stagnationsdrucks
#headers = ['Pstag','G','Paus','ob kritisch']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for P_stag in Pstag_List:
#    Zustand=calc_G_HEM_Zweiphasig(P_stag,T_0,p_2,T_2,zeta)
#    G_List.append(Zustand[0])
#    row=[P_stag]+Zustand
#    rows.append(row)
#with open('Pstag.csv','w',newline ='') as f:  #schreiben in ein Datei
#    f_csv = csv.writer(f)
#    f_csv.writerow(headers)
#    f_csv.writerows(rows)
#plt.figure(1)   #Diagramm zeichnen
#plt.xlabel("Pstag in MPa")
#plt.ylabel("G")
#plt.plot(Pstag_List, G_List)



##unterschiedlicher Umgebungsdruck ohne berücksichtigung des Eitrittsdrucks
#G_List=[]
#Pum_List=np.linspace(0.2,8,21)  #Die Serie des Umgebungsdrucks
#headers = ['Pum','G','Paus','ob kritisch']  #fuer Tabelle
#rows = []  #fuer Tabelle
#for P_um in Pum_List:
#    Zustand=calc_G_HEM_Zweiphasig(p_0,T_0,P_um,T_2,zeta)
#    G_List.append(Zustand[0])
#    row=[P_um]+Zustand
#    rows.append(row)
#with open('Pum.csv','w',newline ='') as g:  #schreiben in ein Datei
#    g_csv = csv.writer(g)
#    g_csv.writerow(headers)
#    g_csv.writerows(rows)
#plt.figure(2)    #Diagramm zeichnen
#plt.xlabel("Pum in MPa")
#plt.ylabel("G")
#plt.plot(Pum_List, G_List)



