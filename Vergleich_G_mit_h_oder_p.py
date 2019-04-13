from iapws import IAPWS97
import math 
import numpy as np

T_0 = 90 + 273.15  #Stagnationstemperatur  K
p_0 = 4           #Stagnationsdruck   MPa
T_2 = 20 + 273.1    #Umgebungstemperatur  K
p_2 = p_0-0.2           #Umgebungsdruck  MPa

state_0 = IAPWS97(T = T_0, P = p_0)
rho_0=state_0.rho
h_0 = state_0.h
s_0 = state_0.s
v_0=state_0.v
state_2= IAPWS97(P=p_2, s=s_0)
h_2=state_2.h
v_2=state_2.v


G_p=math.sqrt(2*(p_0-p_2)*1000000*rho_0)
G_h=math.sqrt(2*(h_0-h_2)*1000/v_2/v_2)

print(p_2,G_p,G_h)
