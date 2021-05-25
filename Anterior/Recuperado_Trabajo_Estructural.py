# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 19:44:33 2020

@author: Santiago Jiménez Gómez
"""

### Trabajo Estructural ###

# Librerias

import sympy as sy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import init_printing
init_printing() 

#%% Funciones 

def radians(x):
    r = (x*sy.pi)/180
    return r

def Matriz_LOC_Portico_P(A,E,L,I): #se ingresa un valor en phi distinto de cero si es el caso de que el empotramiento se encuentre inclinado
    ### El elemento tipo portico plano es la superposición de un elemento cercha plana y uno tipo viga. En dirección del eje X'e 
                    #se comporta como una barra, mientras que en dirección Y'e se comporta como una viga ###

    #sistema de ecuaciones en coordenadas locales
    
    SLoc = np.zeros([6,6])
    
    SLoc[0,0] = (A*E)/L
    SLoc[0,1] = 0
    SLoc[0,2] = 0 
    SLoc[0,3] = (-1*A*E)/L
    SLoc[0,4] = 0 
    SLoc[0,5] = 0
    
    SLoc[1,0] = 0
    SLoc[1,1] = 12*(E*I)/(L**3)
    SLoc[1,2] = 6*(E*I)/(L**2)
    SLoc[1,3] = 0
    SLoc[1,4] = -12*(E*I)/(L**3)
    SLoc[1,5] = 6*(E*I)/(L**2)
    
    SLoc[2,0] = 0 
    SLoc[2,1] = 6*(E*I)/(L**2)
    SLoc[2,2] = 4*(E*I)/(L)
    SLoc[2,3] = 0
    SLoc[2,4] = -6*(E*I)/(L**2)
    SLoc[2,5] = 2*(E*I)/(L)
    
    SLoc[3,0] = (-1*A*E)/L
    SLoc[3,1] = 0
    SLoc[3,2] = 0
    SLoc[3,3] = (A*E)/L
    SLoc[3,4] = 0
    SLoc[3,5] = 0

    SLoc[4,0] = 0
    SLoc[4,1] = -12*(E*I)/(L**3)
    SLoc[4,2] = -6*(E*I)/(L**2)
    SLoc[4,3] = 0
    SLoc[4,4] = 12*(E*I)/(L**3)
    SLoc[4,5] = -6*(E*I)/(L**2)

    SLoc[5,0] = 0
    SLoc[5,1] = 6*(E*I)/(L**2)
    SLoc[5,2] = 2*(E*I)/(L)
    SLoc[5,3] = 0
    SLoc[5,4] = -6*(E*I)/(L**2)
    SLoc[5,5] = 4*(E*I)/(L)
    
    return np.asmatrix(SLoc,dtype = float) 


def Matriz_trans_P_plano(teta,phi):
        #Matriz de transformación 
    
    MTG = np.zeros([6,6])
    
    MTG[0,0] = sy.cos(teta - phi)
    MTG[0,1] = sy.sin(teta - phi)
    MTG[0,2] = 0
    MTG[0,3] = 0
    MTG[0,4] = 0
    MTG[0,5] = 0

    MTG[1,0] = -1*sy.sin(teta - phi)
    MTG[1,1] = sy.cos(teta - phi)
    MTG[1,2] = 0
    MTG[1,3] = 0 
    MTG[1,4] = 0 
    MTG[1,5] = 0

    MTG[2,0] = 0 
    MTG[2,1] = 0 
    MTG[2,2] = 1 
    MTG[2,3] = 0
    MTG[2,4] = 0
    MTG[2,5] = 0

    MTG[3,0] = 0 
    MTG[3,1] = 0 
    MTG[3,2] = 0
    MTG[3,3] = sy.cos(teta - phi)
    MTG[3,4] = sy.sin(teta - phi)
    MTG[3,5] = 0

    MTG[4,0] = 0
    MTG[4,1] = 0
    MTG[4,2] = 0
    MTG[4,3] = -1*sy.sin(teta - phi)
    MTG[4,4] = sy.cos(teta - phi)
    MTG[4,5] = 0

    MTG[5,0] = 0 
    MTG[5,1] = 0
    MTG[5,2] = 0
    MTG[5,3] = 0
    MTG[5,4] = 0
    MTG[5,5] = 1
    
    return np.asmatrix(MTG,dtype = float)

def V_empotramiento_P_plano(p,q,L,x,xi,xj):
        #Funciones de forma 
        
    Psi1=1-(x/L)
    Psi2=1-3*(x/L)**2+2*(x/L)**3
    Psi3=((x/L)-2*(x/L)**2+(x/L)**3)*L
    Psi4=(x/L)
    Psi5=3*(x/L)**2-2*(x/L)**3
    Psi6=(-(x/L)**2+(x/L)**3)*L
    

    
    #Funcion de Green
    
#    Gxx1 = (L/(A*E))*Psi4*Psi1.subs({x:xi})
#    Gxx2 = (L/(A*E))*Psi1*Psi4.subs({x:xi})
#    
#    Gyy1=L**3/(6*E*I)*(-(x/L  )**3*Psi2.subs({x:xi})+3*(x/L  )**2*Psi3.subs({x:xi})/L)
#    Gyy2=L**3/(6*E*I)*(-(1-x/L)**3*Psi5.subs({x:xi})-3*(1-x/L)**2*Psi6.subs({x:xi})/L)
    
    Vemp =np.zeros([6,1])
    
    Vemp[0,0] = -1*sy.integrate(Psi1*p,(x,xi,xj))
    Vemp[1,0] = -1*sy.integrate(Psi2*q,(x,xi,xj))
    Vemp[2,0] = -1*sy.integrate(Psi3*q,(x,xi,xj))
    Vemp[3,0] = -1*sy.integrate(Psi4*p,(x,xi,xj))
    Vemp[4,0] = -1*sy.integrate(Psi5*q,(x,xi,xj))
    Vemp[5,0] = -1*sy.integrate(Psi6*q,(x,xi,xj))
    

    return np.asmatrix(Vemp,dtype = float)


def Matriz_pilas(A,E,L,I,kl):
    
    Lambda=(kl/(4*E*I))**0.25
    
    s =np.sin(Lambda*L)
    c =np.cos(Lambda*L)
    sh=np.sinh(Lambda*L)
    ch=np.cosh(Lambda*L)
    
    #   Funciones de forma
    
    Psi2=(-(s**2*ch**2+c**2*sh**2)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+(s*c+sh*ch)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-(s*c+sh*ch)*sy.cos(Lambda*x)*sy.sinh(Lambda*x)+(sh**2-s**2)*sy.cos(Lambda*x)*sy.cosh(Lambda*x))/(sh**2-s**2)
    Psi3=(1/Lambda)*((s*c-sh*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+sh**2*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s**2*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi5=(2*s*sh*sy.sin(Lambda*x)*sy.sinh(Lambda*x)-(s*ch+c*sh)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)+(s*ch+c*sh)*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi6=(1/Lambda)*((c*sh-s*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+s*sh*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s*sh*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    
    #   Funciones de Green de una viga sobre fundacion flexible doblemente empotrada
    
    Gyy1=1/E*I*(-(sy.sin(Lambda*x)*sy.cosh(Lambda*x)-sy.sinh(Lambda*x)*sy.cos(Lambda*x))/(4*Lambda**3)*Psi2.subs({x:xi})+sy.sin(Lambda*x)*sy.sinh(Lambda*x)/(2*Lambda**2)*Psi3.subs({x:xi}))
    Gyy2=1/E*I*(-(sy.sin(Lambda*(L-x))*sy.cosh(Lambda*(L-x))-sy.sinh(Lambda*(L-x))*sy.cos(Lambda*(L-x)))/(4*Lambda**3)*Psi5.subs({x:xi})-sy.sin(Lambda*(L-x))*sy.sinh(Lambda*(L-x))/(2*Lambda**2)*Psi6.subs({x:xi}))
    
    # Matriz de rigidez
    
    K=np.zeros([6,6])
    
    K[0,0] = A*E/L
    K[0,1] = 0
    K[0,2] = 0
    K[0,3] = -A*E/L
    K[0,4] = 0
    K[0,5] = 0
    
    K[1,0] = 0
    K[1,1] = 4*E*I*Lambda**3*(s*c+sh*ch)/(sh**2-s**2)
    K[1,2] = 2*E*I*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
    K[1,3] = 0
    K[1,4] = -4*E*I*Lambda**3*(s*ch+c*sh)/(sh**2-s**2)
    K[1,5] = 4*E*I*Lambda**2*(s*sh)/(sh**2-s**2)
    
    K[2,0] = 0
    K[2,1] = 2*E*I*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
    K[2,2] = 2*E*I*Lambda*(sh*ch-s*c)/(sh**2-s**2)
    K[2,3] = 0
    K[2,4] = -4*E*I*Lambda**2*(s*sh)/(sh**2-s**2)
    K[2,5] = 2*E*I*Lambda*(s*ch-c*sh)/(sh**2-s**2)
    
    K[3,0] = -A*E/L
    K[3,1] = 0
    K[3,2] = 0
    K[3,3] = A*E/L
    K[3,4] = 0
    K[3,5] = 0

    K[4,0] = 0
    K[4,1] = -4*E*I*Lambda**3*(s*ch+c*sh)/(sh**2-s**2)
    K[4,2] = -4*E*I*Lambda**2*(s*sh)/(sh**2-s**2)
    K[4,3] = 0
    K[4,4] = 4*E*I*Lambda**3*(s*c+sh*ch)/(sh**2-s**2)
    K[4,5] = -2*E*I*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
    
    K[5,0] = 0 
    K[5,1] = 4*E*I*Lambda**2*(s*sh)/(sh**2-s**2) 
    K[5,2] = 2*E*I*Lambda*(s*ch-c*sh)/(sh**2-s**2)
    K[5,3] = 0
    K[5,4] = -2*E*I*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
    K[5,5] = 2*E*I*Lambda*(sh*ch-s*c)/(sh**2-s**2)
    
    
    return np.asmatrix(K,dtype = float)

def Vec_pila_emp(p,q,E,I,kl,L):  #Calculo de las fuerzas de empotramiento
    
    Lambda=(kl/(4*E*I))**0.25
    
    s =np.sin(Lambda*L)
    c =np.cos(Lambda*L)
    sh=np.sinh(Lambda*L)
    ch=np.cosh(Lambda*L)
    
    Psi1=1-(x/L)
    Psi4=(x/L)
    Psi2=(-(s**2*ch**2+c**2*sh**2)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+(s*c+sh*ch)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-(s*c+sh*ch)*sy.cos(Lambda*x)*sy.sinh(Lambda*x)+(sh**2-s**2)*sy.cos(Lambda*x)*sy.cosh(Lambda*x))/(sh**2-s**2)
    Psi3=(1/Lambda)*((s*c-sh*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+sh**2*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s**2*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi5=(2*s*sh*sy.sin(Lambda*x)*sy.sinh(Lambda*x)-(s*ch+c*sh)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)+(s*ch+c*sh)*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi6=(1/Lambda)*((c*sh-s*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+s*sh*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s*sh*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    
    FEmp=np.zeros([6,1])
    
    FEmp[0]=sy.N(-sy.integrate(p*Psi1,(x,0,L)))
    FEmp[1]=sy.N(-sy.integrate(q*Psi2,(x,0,L)))
    FEmp[2]=sy.N(-sy.integrate(q*Psi3,(x,0,L)))
    FEmp[3]=sy.N(-sy.integrate(p*Psi4,(x,0,L)))
    FEmp[4]=sy.N(-sy.integrate(q*Psi5,(x,0,L)))
    FEmp[5]=sy.N(-sy.integrate(q*Psi6,(x,0,L)))
    
    return np.asarray(FEmp,dtype = float)


def campo_desplazamientos_pilas(ui,uj,vi,vj,thetai,thetaj,E,I,L,p,q,A):
    Lambda=(kl/(4*E*I))**0.25
    
    s =np.sin(Lambda*L)
    c =np.cos(Lambda*L)
    sh=np.sinh(Lambda*L)
    ch=np.cosh(Lambda*L)
    
    Psi1=1-(x/L)
    Psi4=(x/L)
    Psi2=(-(s**2*ch**2+c**2*sh**2)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+(s*c+sh*ch)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-(s*c+sh*ch)*sy.cos(Lambda*x)*sy.sinh(Lambda*x)+(sh**2-s**2)*sy.cos(Lambda*x)*sy.cosh(Lambda*x))/(sh**2-s**2)
    Psi3=(1/Lambda)*((s*c-sh*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+sh**2*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s**2*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi5=(2*s*sh*sy.sin(Lambda*x)*sy.sinh(Lambda*x)-(s*ch+c*sh)*sy.sin(Lambda*x)*sy.cosh(Lambda*x)+(s*ch+c*sh)*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    Psi6=(1/Lambda)*((c*sh-s*ch)*sy.sin(Lambda*x)*sy.sinh(Lambda*x)+s*sh*sy.sin(Lambda*x)*sy.cosh(Lambda*x)-s*sh*sy.cos(Lambda*x)*sy.sinh(Lambda*x))/(sh**2-s**2)
    
    #Funcion de Green

    
    
    Gxx1 = (L/(A*E))*Psi4*Psi1.subs({x:xi})
    Gxx2 = (L/(A*E))*Psi1*Psi4.subs({x:xi})
    
    Gyy1=L**3/(6*E*I)*(-(x/L  )**3*Psi2.subs({x:xi})+3*(x/L  )**2*Psi3.subs({x:xi})/L)
    Gyy2=L**3/(6*E*I)*(-(1-x/L)**3*Psi5.subs({x:xi})-3*(1-x/L)**2*Psi6.subs({x:xi})/L)
    
    
    #   Funciones de Green de una viga sobre fundacion flexible doblemente empotrada
    
    Gyy1=1/E*I*(-(sy.sin(Lambda*x)*sy.cosh(Lambda*x)-sy.sinh(Lambda*x)*sy.cos(Lambda*x))/(4*Lambda**3)*Psi2.subs({x:xi})+sy.sin(Lambda*x)*sy.sinh(Lambda*x)/(2*Lambda**2)*Psi3.subs({x:xi}))
    Gyy2=1/E*I*(-(sy.sin(Lambda*(L-x))*sy.cosh(Lambda*(L-x))-sy.sinh(Lambda*(L-x))*sy.cos(Lambda*(L-x)))/(4*Lambda**3)*Psi5.subs({x:xi})-sy.sin(Lambda*(L-x))*sy.sinh(Lambda*(L-x))/(2*Lambda**2)*Psi6.subs({x:xi}))
    
    #Campo Homogeneo
    
    UEh = Psi1*ui + Psi4*uj
    VEh = Psi2*vi + Psi3*thetai + Psi5*vj + Psi6*thetaj
    
    #Campo Empotrado
    
    UEf = sy.integrate(p*Gxx2,(xi,0,x)) + sy.integrate(p*Gxx1,(xi,x,L))
    VEf = sy.integrate(q*Gyy2,(xi,0,x)) + sy.integrate(q*Gyy1,(xi,x,L))
    
    U_E = UEh + UEf
    V_E = VEh + VEf

    return U_E,V_E
    
def campo_desplazamientos_portico(ui,uj,vi,vj,thetai,thetaj,E,I,L,p,q,A):
    
    Psi1=1-(x/L)
    Psi2=1-3*(x/L)**2+2*(x/L)**3
    Psi3=((x/L)-2*(x/L)**2+(x/L)**3)*L
    Psi4=(x/L)
    Psi5=3*(x/L)**2-2*(x/L)**3
    Psi6=(-(x/L)**2+(x/L)**3)*L    
    
    #Funcion de Green

    Gxx1 = (L/(A*E))*Psi4*Psi1.subs({x:xi})
    Gxx2 = (L/(A*E))*Psi1*Psi4.subs({x:xi})
    
    Gyy1=(L**3/(6*E*I))*(-(x/L  )**3*Psi2.subs({x:xi})+3*(x/L  )**2*Psi3.subs({x:xi})/L)
    Gyy2=(L**3/(6*E*I))*(-(1-(x/L))**3*Psi5.subs({x:xi})-3*(1-(x/L))**2*Psi6.subs({x:xi})/L)
    
    
    #Campo Homogeneo
    
    UEh = Psi1*ui + Psi4*uj
    VEh = Psi2*vi + Psi3*thetai + Psi5*vj + Psi6*thetaj
    
    #Campo Empotrado
    
    UEf = sy.integrate(p.subs({x:xi})*Gxx2,(xi,0,x)) + sy.integrate(p.subs({x:xi})*Gxx1,(xi,x,L))
    VEf = sy.integrate(q.subs({x:xi})*Gyy2,(xi,0,x)) + sy.integrate(q.subs({x:xi})*Gyy1,(xi,x,L))
    
    
    U_E = UEh + UEf
    V_E = VEh + VEf
    
    return U_E,V_E
    
#%% 1. Discretizar Elemetos - Nodos - Ejes

x = sy.symbols('x')
xi = sy.symbols('xi')

a = 3
B = 5
c = 3
d = 3
e = 3

kl = 5100

Q1 = 42
Q2 = 35 

E = 19000000 

    # Calculo de intercia para elementos 
    # pilas
    
r = 0.45
A_PI = np.pi*(r**2)
Ix_pi = (1/4)*np.pi*(r**4) #Inercia para elemento tipo pila

    #Porticos planos de sección rectangular
    
b = 0.35
h = 0.40     
A_P = b*h
Ix_P = (1/12)*b*(h**3) #Inercia para elementos tipo portico

#%% Elementos tipo Pila

#   Elemento A

LA = d
Theta_A = radians(90)
Phi_A  = 0
Mloc_A = Matriz_pilas(A_PI,E,LA,Ix_pi,kl)
Mtrans_A = (Matriz_trans_P_plano(Theta_A,Phi_A))
Vemp_A = Vec_pila_emp(0,0,E,Ix_pi,kl,LA)

Mtranspose_A = np.transpose(Mtrans_A)
KA = Mtranspose_A @ Mloc_A @ Mtrans_A
Vemp_GA = np.transpose(Mtrans_A)@Vemp_A

#   Elemento E

LE = d
Theta_E = radians(90)
Phi_E  = 0
Mloc_E = Matriz_pilas(A_PI,E,LE,Ix_pi,kl)
Mtrans_E = (Matriz_trans_P_plano(Theta_E,Phi_E))
Vemp_E = Vec_pila_emp(0,0,E,Ix_pi,kl,LE)

Mtranspose_E = np.transpose(Mtrans_E)
KE = Mtranspose_E @ Mloc_E @ Mtrans_E
Vemp_GE = np.transpose(Mtrans_E)@Vemp_E


#%% Elementos tipo Portico Plano

#   Elemento B

LB = sy.sqrt(a**2 + e**2)
Theta_B = sy.atan2(3,3)
Phi_B  = 0

QB_Globales = (-Q1/LB)*x  ## Hacia Arriba
p_B = QB_Globales*sy.Abs(sy.cos(Theta_B))*sy.sin(Theta_B)
q_B = QB_Globales*sy.Abs(sy.cos(Theta_B))*sy.cos(Theta_B)


Mloc_B = (Matriz_LOC_Portico_P(A_P,E,LB,Ix_P))
Mtrans_B = (Matriz_trans_P_plano(Theta_B,Phi_B))
Vemp_B = V_empotramiento_P_plano(p_B,q_B,LB,x,0,LB)


Mtranspose_B = np.transpose(Mtrans_B)
KB = Mtranspose_B @ Mloc_B @ Mtrans_B
Vemp_GB = np.transpose(Mtrans_B)@Vemp_B

#   Elemento C

LC = B
Theta_C = 0
Phi_C  = 0

QC_Globales = (Q1/(LC/2))*x - Q1 ## Hacia Arriba
p_C = 0*x
q_C1 = QC_Globales
q_C2 = 0*x

Mloc_C = (Matriz_LOC_Portico_P(A_P,E,LC,Ix_P))
Mtrans_C = Matriz_trans_P_plano(Theta_C,Phi_C)
Vemp_C = V_empotramiento_P_plano(p_C,q_C1,LC,x,0,LC/2)

Mtranspose_C = (np.transpose(Mtrans_C))
KC = Mtranspose_C @ Mloc_C @ Mtrans_C
Vemp_GC = np.transpose(Mtrans_C)@Vemp_C

#   Elemento D

LD = sy.sqrt(e**2 + c**2)
Theta_D = -sy.atan2(3,3)
Phi_D  = 0

QD_Globales = (Q2/c)*x - Q2 ## Hacia Afuera
p_D = QD_Globales*sy.Abs(sy.sin(Theta_D))*sy.cos(Theta_D)
q_D = -QD_Globales*sy.Abs(sy.sin(Theta_D))*sy.sin(Theta_D)


Mloc_D = (Matriz_LOC_Portico_P(A_P,E,LD,Ix_P))
Mtrans_D = (Matriz_trans_P_plano(Theta_D,Phi_D))
Vemp_D = V_empotramiento_P_plano(p_D,q_D,LD,x,0,LD)

Mtranspose_D = (np.transpose(Mtrans_D))
KD = Mtranspose_D @ Mloc_D @ Mtrans_D
Vemp_GD = np.transpose(Mtrans_D)@Vemp_D


#%% Sistema matricial

DesNod = np.zeros([15,15])

# Fila 1 ### M1
DesNod[0,0] = KA[2,2]
DesNod[0,1] = KA[2,3]
DesNod[0,2] = KA[2,4]
DesNod[0,3] = KA[2,5]
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] = 
#DesNod[0,7] = 
#DesNod[0,8] = 
#DesNod[0,9] = 
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 2 ### FX2
DesNod[1,0] = KA[3,2]
DesNod[1,1] = KA[3,3] + KB[0,0]
DesNod[1,2] = KA[3,4] + KB[0,1]
DesNod[1,3] = KA[3,5] + KB[0,2]
DesNod[1,4] = KB[0,3]
DesNod[1,5] = KB[0,4] 
DesNod[1,6] = KB[0,5]
#DesNod[0,7] = 
#DesNod[0,8] = 
#DesNod[0,9] = 
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 3 ### FY2

DesNod[2,0] = KA[4,2]
DesNod[2,1] = KA[4,3] + KB[1,0]
DesNod[2,2] = KA[4,4] + KB[1,1]
DesNod[2,3] = KA[4,5] + KB[1,2]
DesNod[2,4] = KB[1,3]
DesNod[2,5] = KB[1,4]
DesNod[2,6] = KB[1,5]
#DesNod[1,7] = 
#DesNod[1,8] = 
#DesNod[1,9] = 
#DesNod[1,10] = 
#DesNod[1,11] = 
#DesNod[1,12] =
#DesNod[1,13] = 
#DesNod[1,14] = 

# Fila 4 ### M2

DesNod[3,0] = KA[5,2]               #theta1
DesNod[3,1] = KA[5,3] + KB[2,0]     #u2
DesNod[3,2] = KA[5,4] + KB[2,1]     #v2
DesNod[3,3] = KA[5,5] + KB[2,2]     #theta2
DesNod[3,4] = KB[2,3]               #u3
DesNod[3,5] = KB[2,4]               #v3
DesNod[3,6] = KB[2,5]               #theta3
#DesNod[0,7] = 
#DesNod[0,8] = 
#DesNod[0,9] = 
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 5 ### FX3

#DesNod[0,0] = 
DesNod[4,1] = KB[3,0]
DesNod[4,2] = KB[3,1]
DesNod[4,3] = KB[3,2]
DesNod[4,4] = KB[3,3] + KC[0,0]
DesNod[4,5] = KB[3,4] + KC[0,1]
DesNod[4,6] = KB[3,5] + KC[0,2]
DesNod[4,7] = KC[0,3]
DesNod[4,8] = KC[0,4]
DesNod[4,9] = KC[0,5]
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 6 ### FY3
#DesNod[0,0] = 
DesNod[5,1] = KB[4,0]               #u2
DesNod[5,2] = KB[4,1]               #v2
DesNod[5,3] = KB[4,2]               #theta2
DesNod[5,4] = KB[4,3] + KC[1,0]     #u3
DesNod[5,5] = KB[4,4] + KC[1,1]     #v3
DesNod[5,6] = KB[4,5] + KC[1,2]     #theta3
DesNod[5,7] = KC[1,3]               #u4
DesNod[5,8] = KC[1,4]               #v4
DesNod[5,9] = KC[1,5]               #theta4C
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 7 ### M3
#DesNod[0,0] = 
DesNod[6,1] = KB[5,0]           #u2
DesNod[6,2] = KB[5,1]           #v2
DesNod[6,3] = KB[5,2]           #theta2
DesNod[6,4] = KB[5,3] + KC[2,0] #u3
DesNod[6,5] = KB[5,4] + KC[2,1] #v3
DesNod[6,6] = KB[5,5] + KC[2,2] #theta3
DesNod[6,7] = KC[2,3]           #u4
DesNod[6,8] = KC[2,4]           #v4
DesNod[6,9] = KC[2,5]           #theta4C
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 8 ### FX4

#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
DesNod[7,4] = KC[3,0]               #u3
DesNod[7,5] = KC[3,1]               #v3
DesNod[7,6] = KC[3,2]               #theta3
DesNod[7,7] = KC[3,3] + KD[0,0]     #u4
DesNod[7,8] = KC[3,4] + KD[0,1]     #v4
DesNod[7,9] = KC[3,5]               # theta 4C
DesNod[7,10] = KD[0,2]              # theta 4D
DesNod[7,11] = KD[0,3]              #u5
DesNod[7,12] = KD[0,4]              #v5
DesNod[7,13] = KD[0,5]              #theta5
#DesNod[0,14] = 

# Fila 9 ### FY4

#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] =                
DesNod[8,4] = KC[4,0]              #u3
DesNod[8,5] = KC[4,1]               #v3
DesNod[8,6] = KC[4,2]               #theta3
DesNod[8,7] = KC[4,3] + KD[1,0]     #u4
DesNod[8,8] = KC[4,4] + KD[1,1]     #v4
DesNod[8,9] = KC[4,5]               # theta 4C
DesNod[8,10] = KD[1,2]              # theta 4D
DesNod[8,11] = KD[1,3]              #u5
DesNod[8,12] = KD[1,4]              #v5
DesNod[8,13] = KD[1,5]              #theta5
#DesNod[0,14] = 

# Fila 10 ### M4C
#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
DesNod[9,4] = KC[5,0]
DesNod[9,5] = KC[5,1]
DesNod[9,6] = KC[5,2]
DesNod[9,7] = KC[5,3]
DesNod[9,8] = KC[5,4]
DesNod[9,9] = KC[5,5]
#DesNod[0,10] = 
#DesNod[0,11] = 
#DesNod[0,12] =
#DesNod[0,13] = 
#DesNod[0,14] = 

# Fila 11 ### M4D
#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] =
DesNod[10,7] = KD[2,0]
DesNod[10,8] = KD[2,1]
DesNod[10,9] = 0            #theta 4C
DesNod[10,10] = KD[2,2]      #theta 4D 
DesNod[10,11] = KD[2,3]
DesNod[10,12] = KD[2,4]
DesNod[10,13] = KD[2,5]
#DesNod[0,14] = 

# Fila 12 ### FX5

#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] = 
DesNod[11,7] = KD[3,0]              #u4
DesNod[11,8] = KD[3,1]              #v4
DesNod[11,9] = 0                    #theta4C
DesNod[11,10] = KD[3,2]             #theta4D
DesNod[11,11] = KD[3,3] + KE[3,3]   #u5
DesNod[11,12] = KD[3,4] + KE[3,4]   #v5
DesNod[11,13] = KD[3,5] + KE[3,5]   #theta5
DesNod[11,14] = KE[3,2]             #theta6

# Fila 13 ### FY5

#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] = 
DesNod[12,7] = KD[4,0]
DesNod[12,8] = KD[4,1]
DesNod[12,9] = 0
DesNod[12,10] = KD[4,2]
DesNod[12,11] = KD[4,3] + KE[4,3]
DesNod[12,12] = KD[4,4] + KE[4,4]
DesNod[12,13] = KD[4,5] + KE[4,5]
DesNod[12,14] = KE[4,2]

# Fila 14 ### M5

#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] = 
DesNod[13,7] = KD[5,0]
DesNod[13,8] = KD[5,1]
DesNod[13,9] = 0
DesNod[13,10] = KD[5,2]
DesNod[13,11] = KD[5,3] + KE[5,3]
DesNod[13,12] = KD[5,4] + KE[5,4]
DesNod[13,13] = KD[5,5] + KE[5,5]
DesNod[13,14] = KE[5,2]

# Fila 15 ### M6
#DesNod[0,0] = 
#DesNod[0,1] = 
#DesNod[0,2] = 
#DesNod[0,3] = 
#DesNod[0,4] = 
#DesNod[0,5] = 
#DesNod[0,6] = 
#DesNod[0,7] = 
#DesNod[0,8] = 
#DesNod[0,9] = 
#DesNod[0,10] = 
DesNod[14,11] = KE[2,3] 
DesNod[14,12] = KE[2,4]
DesNod[14,13] = KE[2,5]
DesNod[14,14] = KE[2,2]

FUE_EMP = np.zeros([15,1])

FUE_EMP[0,0] = Vemp_GA[2,0]                      #theta1
FUE_EMP[1,0] = Vemp_GA[3,0] + Vemp_GB[0,0]        #u2
FUE_EMP[2,0] = Vemp_GA[4,0] + Vemp_GB[1,0]        #v2
FUE_EMP[3,0] = Vemp_GA[5,0] + Vemp_GB[2,0]        #theta2
FUE_EMP[4,0] = Vemp_GB[3,0] + Vemp_GC[0,0]        #u3
FUE_EMP[5,0] = Vemp_GB[4,0] + Vemp_GC[1,0]        #v3
FUE_EMP[6,0] = Vemp_GB[5,0] + Vemp_GC[2,0]        #theta3
FUE_EMP[7,0] = Vemp_GC[3,0] + Vemp_GD[0,0]        #u4
FUE_EMP[8,0] = Vemp_GC[4,0] + Vemp_GD[1,0]        #v4
FUE_EMP[9,0] = Vemp_GC[5,0]                      #theta4C
FUE_EMP[10,0] = Vemp_GD[2,0]                     #theta4D
FUE_EMP[11,0] = Vemp_GD[3,0] + Vemp_GE[3,0]       #u5
FUE_EMP[12,0] = Vemp_GD[4,0] + Vemp_GE[4,0]       #v5
FUE_EMP[13,0] = Vemp_GD[5,0] + Vemp_GE[5,0]       #theta5
FUE_EMP[14,0] = Vemp_GE[2,0]                     #theta6
 
Desplazamientos = np.linalg.solve(DesNod,(-1*FUE_EMP))


#%% notacion

######## Desplazamientos en notacion cientifica ####################
l = []
for i in range(15):
    v = np.format_float_scientific(Desplazamientos[i,0], precision=2,exp_digits=1)
    l.append(v)
    
des_notacion = (sy.Matrix(np.transpose(np.asmatrix(l))))

####################################################################

u1 = 0
v1 = 0
theta1 = Desplazamientos[0,0]
u2 = Desplazamientos[1,0]
v2 = Desplazamientos[2,0]
theta2 = Desplazamientos[3,0]
u3 = Desplazamientos[4,0]
v3 = Desplazamientos[5,0]
theta3 = Desplazamientos[6,0]
u4 = Desplazamientos[7,0]
v4 = Desplazamientos[8,0]
theta4C = Desplazamientos[9,0]
theta4D = Desplazamientos[10,0]
u5 = Desplazamientos[11,0]
v5 = Desplazamientos[12,0]
theta5 = Desplazamientos[13,0]
theta6 = Desplazamientos[14,0]
u6 = 0
v6 = 0


#%%     Cálculo de los desplazamientos nodales en coordenadas locales

#############################################################

#       Elemento A

DesNod_Glo_A = np.zeros([6,1])

DesNod_Glo_A[0,0]=u1
DesNod_Glo_A[1,0]=v1
DesNod_Glo_A[2,0]=theta1
DesNod_Glo_A[3,0]=u2
DesNod_Glo_A[4,0]=v2
DesNod_Glo_A[5,0]=theta2

DesNod_Loc_A = Mtrans_A@DesNod_Glo_A

u1_A = DesNod_Loc_A[0,0]
v1_A = DesNod_Loc_A[1,0]
theta1_A = DesNod_Loc_A[2,0]
u2_A = DesNod_Loc_A[3,0]
v2_A = DesNod_Loc_A[4,0]
theta2_A = DesNod_Loc_A[5,0]

#       Elemento E

DesNod_Glo_E=np.zeros([6,1])

DesNod_Glo_E[0,0]=u6
DesNod_Glo_E[1,0]=v6
DesNod_Glo_E[2,0]=theta6
DesNod_Glo_E[3,0]=u5
DesNod_Glo_E[4,0]=v5
DesNod_Glo_E[5,0]=theta5

DesNod_Loc_E = Mtrans_E@DesNod_Glo_E

u6_E = DesNod_Loc_E[0,0]
v6_E = DesNod_Loc_E[1,0]
theta6_E = DesNod_Loc_E[2,0]
u5_E = DesNod_Loc_E[3,0]
v5_E = DesNod_Loc_E[4,0]
theta5_E = DesNod_Loc_E[5,0]


################################################################

#       Elemento B

DesNod_Glo_B = np.zeros([6,1])

DesNod_Glo_B[0,0]=u2
DesNod_Glo_B[1,0]=v2
DesNod_Glo_B[2,0]=theta2
DesNod_Glo_B[3,0]=u3
DesNod_Glo_B[4,0]=v3
DesNod_Glo_B[5,0]=theta3

DesNod_Loc_B = Mtrans_B@DesNod_Glo_B

u2_B = DesNod_Loc_B[0,0]
v2_B = DesNod_Loc_B[1,0]
theta2_B = DesNod_Loc_B[2,0]
u3_B = DesNod_Loc_B[3,0]
v3_B = DesNod_Loc_B[4,0]
theta3_B = DesNod_Loc_B[5,0]

#       Elemento C

DesNod_Glo_C=np.zeros([6,1])

DesNod_Glo_C[0,0]=u3
DesNod_Glo_C[1,0]=v3
DesNod_Glo_C[2,0]=theta3
DesNod_Glo_C[3,0]=u4
DesNod_Glo_C[4,0]=v4
DesNod_Glo_C[5,0]=theta4C

DesNod_Loc_C = Mtrans_C@DesNod_Glo_C

u3_C = DesNod_Loc_C[0,0]
v3_C = DesNod_Loc_C[1,0]
theta3_C = DesNod_Loc_C[2,0]
u4_C = DesNod_Loc_C[3,0]
v4_C = DesNod_Loc_C[4,0]
theta4_C = DesNod_Loc_C[5,0]

#       Elemento D

DesNod_Glo_D=np.zeros([6,1])

DesNod_Glo_D[0,0]=u4
DesNod_Glo_D[1,0]=v4
DesNod_Glo_D[2,0]=theta4D
DesNod_Glo_D[3,0]=u5
DesNod_Glo_D[4,0]=v5
DesNod_Glo_D[5,0]=theta5

DesNod_Loc_D = Mtrans_D@DesNod_Glo_D

u4_D = DesNod_Loc_D[0,0]
v4_D = DesNod_Loc_D[1,0]
theta4_D = DesNod_Loc_D[2,0]
u5_D = DesNod_Loc_D[3,0]
v5_D = DesNod_Loc_D[4,0]
theta5_D = DesNod_Loc_D[5,0]

#%% Calculo de las Reacciones

Vector_Des_A = np.zeros([6,1])
Vector_Des_A[0,0] = 0
Vector_Des_A[1,0] = 0
Vector_Des_A[2,0] = theta1
Vector_Des_A[3,0] = u2
Vector_Des_A[4,0] = v2
Vector_Des_A[5,0] = theta2

Vector_Des_E = np.zeros([6,1])
Vector_Des_E[0,0] = 0
Vector_Des_E[1,0] = 0
Vector_Des_E[2,0] = theta6
Vector_Des_E[3,0] = u5
Vector_Des_E[4,0] = v5
Vector_Des_E[5,0] = theta5

FX1 = KA[0]@Vector_Des_A + Vemp_GA[0,0]
FY1 = KA[1]@Vector_Des_A + Vemp_GA[1,0]

FX6 = KA[0]@Vector_Des_E + Vemp_GE[0,0]
FY6 = KA[1]@Vector_Des_E + Vemp_GE[1,0]


#%% Campos de Desplazamientos en locales

U_A,V_A = campo_desplazamientos_pilas(u1_A,u2_A,v1_A,v2_A,theta1_A,theta2_A,E,Ix_pi,LA,0,0,A_PI)
U_E,V_E = campo_desplazamientos_pilas(u6_E,u5_E,v6_E,v5_E,theta6_E,theta5_E,E,Ix_pi,LE,0,0,A_PI)    
                                        
U_B,V_B = campo_desplazamientos_portico(u2_B,u3_B,v2_B,v3_B,theta2_B,theta3_B,E,Ix_P,LB,p_B,q_B,A_P)
U_C,V_C_BAD = campo_desplazamientos_portico(u3_C,u4_C,v3_C,v4_C,theta3_C,theta4_C,E,Ix_P,LC,p_C,q_C1,A_P)
U_D,V_D = campo_desplazamientos_portico(u4_D,u5_D,v4_D,v5_D,theta4_D,theta5_D,E,Ix_P,LD,p_D,q_D,A_P)

# Elemento C

#q1 va hasta LC/2

Psi1=1-(x/LC)
Psi2=1-3*(x/LC)**2+2*(x/LC)**3
Psi3=((x/LC)-2*(x/LC)**2+(x/LC)**3)*LC
Psi4=(x/LC)
Psi5=3*(x/LC)**2-2*(x/LC)**3
Psi6=(-(x/LC)**2+(x/LC)**3)*LC

#Funcion de Green

Gxx1 = (LC/(A_P*E))*Psi4*Psi1.subs({x:xi})
Gxx2 = (LC/(A_P*E))*Psi1*Psi4.subs({x:xi})

Gyy1=(LC**3/(6*E*Ix_P))*(-(x/LC  )**3*Psi2.subs({x:xi})+3*(x/LC  )**2*Psi3.subs({x:xi})/LC)
Gyy2=(LC**3/(6*E*Ix_P))*(-(1-(x/LC))**3*Psi5.subs({x:xi})-3*(1-(x/LC))**2*Psi6.subs({x:xi})/LC)
    

#Campo Homogeneo

VCh = Psi2*v3_C + Psi3*theta3_C + Psi5*v4_C + Psi6*theta4_C

#Campo Empotrado

VCf1 = sy.integrate(q_C1.subs({x:xi})*Gyy2,(xi,0,x)) + sy.integrate(q_C1.subs({x:xi})*Gyy1,(xi,x,LC/2))
VCf2 = sy.integrate(q_C1.subs({x:xi})*Gyy2,(xi,0,LC/2)) + sy.integrate((0*x).subs({x:xi})*Gyy2,(xi,LC/2,x)) + sy.integrate((0*x).subs({x:xi})*Gyy1,(xi,x,LC))

# Campo Total

V_C1 = VCh + VCf1
V_C2 = VCh + VCf2

# Cálculo de la fuerza que el suelo ejerce sobre las pilas

# Elemento A

F_soil_A=-kl*V_A

# Elemento E

F_soil_E=-kl*V_E

#%% Cálculo de las fuerzas internas 

#       Elemento A
PA = A_PI*E*sy.diff(U_A,x,1)
MA = E*Ix_pi*sy.diff(V_A,x,2)
VA = -E*Ix_pi*sy.diff(V_A,x,3)

#       Elemento E
PE = A_PI*E*sy.diff(U_E,x,1)
ME = E*Ix_pi*sy.diff(V_E,x,2)
VE = -E*Ix_pi*sy.diff(V_E,x,3)


#       Elemento B
PB = A_P*E*sy.diff(U_B,x,1)
MB = E*Ix_P*sy.diff(V_B,x,2)
VB = -E*Ix_P*sy.diff(V_B,x,3)

#       Elemento C
PC = A_P*E*sy.diff(U_C,x,1)

MC1 = E*Ix_P*sy.diff(V_C1,x,2) #        Para x entre 0 y LC/2 
MC2 = E*Ix_P*sy.diff(V_C2,x,2) #        Para x entre LC/2 y LC

VC1 = -E*Ix_P*sy.diff(V_C1,x,3) #       Para x entre 0 y LC/2
VC2 = -E*Ix_P*sy.diff(V_C2,x,3) #       Para x entre LC/2 y LC

#       Elemento D
PD = A_P*E*sy.diff(U_D,x,1)
MD = E*Ix_P*sy.diff(V_D,x,2)
VD = -E*Ix_P*sy.diff(V_D,x,3)


#%% Revision del equilibrio de la estructura



SFX=FX1[0,0]+FX6[0,0]-(Q2*e/2)-sy.integrate(F_soil_A,(x,0,LA))+sy.integrate(F_soil_E,(x,0,LE))
SFY=FY1[0,0]+FY6[0,0]-(Q1*a/2) - (Q1*B/2)/2
#SM1=M1+M3+FY3*8-30*2*2-40*2/2*(3-2/3)+50*(4+2/3)-50*(8-2/3)

#%% Revisiones
#   13.0 Revisiones 

Error=sy.zeros(50,1)

#   13.1 Cumplimiento de ecuaciones diferenciales gobernantes de cada elemento

 #Cumplimiento ecuacion diferencial elemento A
Error[0,0]=sy.expand(A_PI*E*sy.diff(U_A,x,2)+0)
Error[1,0]=sy.expand(Ix_pi*E*sy.diff(V_A,x,4)+(kl*V_A)-0)
 
#Cumplimiento ecuacion diferencial elemento B

Error[2,0]=sy.expand(A_P*E*sy.diff(U_B,x,2)+p_B)
Error[3,0]=sy.expand(Ix_P*E*sy.diff(V_B,x,4)-q_B)

#Cumplimiento ecuacion diferencial elemento C primer tramo 
Error[4,0]=sy.expand(A_P*E*sy.diff(U_C,x,2)+0)
Error[5,0]=sy.expand(Ix_P*E*sy.diff(V_C1,x,4)-q_C1)

#Cumplimiento ecuacion diferencial elemento C segundo tramo
Error[6,0]=sy.expand(A_P*E*sy.diff(U_C,x,2)+0)
Error[7,0]=sy.expand(Ix_P*E*sy.diff(V_C2,x,4)-q_C2)

 #Cumplimiento ecuacion diferencial elemento D
Error[8,0]=sy.expand(A_P*E*sy.diff(U_D,x,2)+p_D)
Error[9,0]=sy.expand(Ix_P*E*sy.diff(V_D,x,4)-q_D)
 
 #Cumplimiento ecuacion diferencial elemento E
Error[10,0]=sy.expand(A_PI*E*sy.diff(U_E,x,2)+0)
Error[11,0]=sy.expand(Ix_pi*E*sy.diff(V_E,x,4)+(kl*V_E)-0)

 
#   13.2 Cumplimiento de continuidad de desplazamientos

#   Nudo 1
Error[12,0]=U_A.subs({x:0})                                                                                                                   
Error[13,0]=V_A.subs({x:0})


#   Nudo 2
Error[14,0]=sy.expand(U_A.subs({x:LA})-U_B.subs({x:0}))
Error[15,0]=sy.expand(V_A.subs({x:LA})-V_B.subs({x:0}))
Error[16,0]=sy.expand(sy.diff(PA,x,1).subs({x:LA})-sy.diff(PB,x,1).subs({x:0}))  
Error[17,0]=sy.expand(sy.diff(VA,x,2).subs({x:LA})-sy.diff(VB,x,2).subs({x:0}))  
Error[18,0]=sy.expand(-sy.diff(MA,x,3).subs({x:LA})+sy.diff(MB,x,3).subs({x:0}))  
 
#   Nudo 3
Error[19,0]=sy.expand(U_B.subs({x:LB})-U_C.subs({x:0}))
Error[20,0]=sy.expand(V_B.subs({x:LB})-V_C1.subs({x:0}))
Error[21,0]=sy.expand(sy.diff(PB,x,1).subs({x:LB})-sy.diff(PC,x,1).subs({x:0}))  
Error[22,0]=sy.expand(sy.diff(VB,x,2).subs({x:LB})-sy.diff(VC1,x,2).subs({x:0}))  
Error[23,0]=sy.expand(-sy.diff(MB,x,3).subs({x:LB})+sy.diff(MC1,x,3).subs({x:0}))   
 
#   Nudo 4
Error[24,0]=sy.expand(U_C.subs({x:LC})-U_D.subs({x:0}))
Error[25,0]=sy.expand(V_C2.subs({x:LC})-V_D.subs({x:0}))
Error[26,0]=sy.expand(sy.diff(PC,x,1).subs({x:LC})-sy.diff(PD,x,1).subs({x:0}))  
Error[27,0]=sy.expand(sy.diff(VC2,x,2).subs({x:LC})-sy.diff(VD,x,2).subs({x:0}))  
Error[28,0]=sy.expand(-sy.diff(MC2,x,3).subs({x:LC})+sy.diff(MD,x,3).subs({x:0})) 

#   Nudo 5
Error[29,0]=sy.expand(U_D.subs({x:LD})-U_E.subs({x:0}))
Error[30,0]=sy.expand(V_D.subs({x:LD})-V_E.subs({x:0}))
Error[31,0]=sy.expand(sy.diff(PD,x,1).subs({x:LD})-sy.diff(PE,x,1).subs({x:0}))  
Error[32,0]=sy.expand(sy.diff(VD,x,2).subs({x:LD})-sy.diff(VE,x,2).subs({x:0}))  
Error[33,0]=sy.expand(-sy.diff(MD,x,3).subs({x:LD})+sy.diff(ME,x,3).subs({x:0}))    
 
#   Nudo 1
Error[34,0]=U_E.subs({x:0})                                                                                                                   
Error[35,0]=V_E.subs({x:0}) 



if SFX or SFY != 0:
    print("NO hay EQUILIBRIO en la estructura, sfx, sfy distintos de 0 \n" "SFX=",  SFX, "SFY=",SFY )
#%% Tablas

#def Tablas_Fuerzas_internas(FI,xi,xj,n,F):
#    import pandas as pd
#    xE = np.linspace(float(xi),float(xj),n)
#    data = {'$x_E$':[], F :[]}
#    for i in xE:
#       valor = np.round(float(FI.subs({x:i})),decimals = 3)
#       data['$x_E$'].append(np.round(float(i),decimals=3))
#       data[str(F)].append(valor)
#    tabla = pd.DataFrame.from_dict(data)
#    
#    return tabla 
#    
#PA_Tabla = Tablas_Fuerzas_internas(PA,0,LA,10,'P(x)')
#VA_Tabla = Tablas_Fuerzas_internas(VA,0,LA,10,'V(x)')
#MA_Tabla = Tablas_Fuerzas_internas(MA,0,LA,10,'M(x)')
#
#PB_Tabla = Tablas_Fuerzas_internas(PB,0,LB,10,'P(x)')
#VB_Tabla = Tablas_Fuerzas_internas(VB,0,LB,10,'V(x)')
#MB_Tabla = Tablas_Fuerzas_internas(MB,0,LB,10,'M(x)')
#
#PD_Tabla = Tablas_Fuerzas_internas(PD,0,LD,10,'P(x)')
#VD_Tabla = Tablas_Fuerzas_internas(VD,0,LD,10,'V(x)')
#MD_Tabla = Tablas_Fuerzas_internas(MD,0,LD,10,'M(x)')
#
#PE_Tabla = Tablas_Fuerzas_internas(PE,0,LE,10,'P(x)')
#VE_Tabla = Tablas_Fuerzas_internas(VE,0,LE,10,'V(x)')
#ME_Tabla = Tablas_Fuerzas_internas(ME,0,LE,10,'M(x)')
#
#PC_Tabla = Tablas_Fuerzas_internas(PC,0,LC,10,'P(x)')
#
#VC1_Tabla = Tablas_Fuerzas_internas(VC1,0,(LC),10,'V(x)').drop(index=[5,6,7,8,9])
#VC2_Tabla = Tablas_Fuerzas_internas(VC2,0,LC,10,'V(x)').drop(index=[0,1,2,3,4])
#
#MC1_Tabla = Tablas_Fuerzas_internas(MC1,0,(LC),10,'M(x)').drop(index=[5,6,7,8,9])
#MC2_Tabla = Tablas_Fuerzas_internas(MC2,0,LC,10,'M(x)').drop(index=[0,1,2,3,4])
#
#tablas = [PA_Tabla,VA_Tabla,MA_Tabla,
#          PB_Tabla,VB_Tabla,MB_Tabla,
#          PD_Tabla,VD_Tabla,MD_Tabla,
#          PE_Tabla,VE_Tabla,ME_Tabla,
#          PC_Tabla,VC1_Tabla,VC2_Tabla,MC1_Tabla,MC2_Tabla]
#
#VC_Tabla = pd.concat([VC1_Tabla,VC2_Tabla], axis=0, sort=False).reset_index(drop=True)
#MC_Tabla = pd.concat([MC1_Tabla,MC2_Tabla], axis=0, sort=False).reset_index(drop=True)
#
#Elemento_A = pd.merge(PA_Tabla,pd.merge(VA_Tabla,MA_Tabla,on='$x_E$'), on='$x_E$')
#Elemento_B = pd.merge(PB_Tabla,pd.merge(VB_Tabla,MB_Tabla,on='$x_E$'), on='$x_E$')
#Elemento_C = pd.merge(PC_Tabla,pd.merge(VC_Tabla,MC_Tabla,on='$x_E$'), on='$x_E$')
#Elemento_D = pd.merge(PD_Tabla,pd.merge(VD_Tabla,MD_Tabla,on='$x_E$'), on='$x_E$')
#Elemento_E = pd.merge(PE_Tabla,pd.merge(VE_Tabla,ME_Tabla,on='$x_E$'), on='$x_E$')
#
#
#Elemento_A.name = 'Elemento_A'
#Elemento_B.name = 'Elemento_B'
#Elemento_C.name = 'Elemento_C'
#Elemento_D.name = 'Elemento_D'
#Elemento_E.name = 'Elemento_E'
#
#Elementos=[Elemento_A,Elemento_B,Elemento_C,Elemento_D,Elemento_E]
#
##for j in tablas:
##    print(j.keys(), '\n' ,j.to_latex(index=False) )
##    
#for jj in Elementos:
#    print(jj.name, '\n' ,jj.to_latex(index=False) )  
    
#%% Figuras

plt.close('all')

EscalaP=0.009
EscalaV=0.03
EscalaM=0.0099

Escala=100
Nx=100

x1,y1 = 0,0
x2,y2 = 0,3
x3,y3 = 3,6
x4,y4 = 8,6
x5,y5 = 11,3
x6,y6 = 11,0

xANum=np.linspace(0,float(LA),Nx)
xBNum=np.linspace(0,float(LB),Nx)
xCNum=np.linspace(0,float(LC),Nx)
#xC2Num=np.linspace(LC/2,float(LC),Nx)
xDNum=np.linspace(0,float(LD),Nx)
xENum=np.linspace(0,float(LE),Nx)

uANum=np.zeros([Nx])
uBNum=np.zeros([Nx])
uCNum=np.zeros([Nx])
#uC2Num=np.zeros([Nx])
uDNum=np.zeros([Nx])
uENum=np.zeros([Nx])

vANum=np.zeros([Nx])
vBNum=np.zeros([Nx])
vCNum=np.zeros([Nx])
#vC2Num=np.zeros([Nx])
vDNum=np.zeros([Nx])
vENum=np.zeros([Nx])

#######################     Campos de fuerzas internas

PANum=np.zeros([Nx])
VANum=np.zeros([Nx])
MANum=np.zeros([Nx])

PENum=np.zeros([Nx])
VENum=np.zeros([Nx])
MENum=np.zeros([Nx])

PBNum=np.zeros([Nx])
VBNum=np.zeros([Nx])
MBNum=np.zeros([Nx])

PCNum=np.zeros([Nx])
VCNum=np.zeros([Nx])#        Para x entre 0 y LC/2 
MCNum=np.zeros([Nx])#        Para x entre 0 y LC/2 
#VC2Num=np.zeros([Nx])#        Para x entre LC/2 y LC
#MC2Num=np.zeros([Nx])#        Para x entre LC/2 y LC

PDNum=np.zeros([Nx])
VDNum=np.zeros([Nx])
MDNum=np.zeros([Nx])




for ix in range(Nx):
    uANum[ix]=U_A.subs({x:xANum[ix]})
    vANum[ix]=V_A.subs({x:xANum[ix]})
    PANum[ix]=PA.subs({x:xANum[ix]})
    VANum[ix]=VA.subs({x:xANum[ix]})
    MANum[ix]=MA.subs({x:xANum[ix]})
    
for ix in range(Nx):
    uBNum[ix]=U_B.subs({x:xBNum[ix]})
    vBNum[ix]=V_B.subs({x:xBNum[ix]})
    PBNum[ix]=PB.subs({x:xBNum[ix]})
    VBNum[ix]=VB.subs({x:xBNum[ix]})
    MBNum[ix]=MB.subs({x:xBNum[ix]})
    
#for ix in range(Nx):
#    uC1Num[ix]=U_C.subs({x:xC1Num[ix]})
#    vC1Num[ix]=V_C1.subs({x:xC1Num[ix]})
#    PCNum[ix] = PC.subs({x:xC1Num[ix]})
#    VC1Num[ix]=VC1.subs({x:xC1Num[ix]})
#    MC1Num[ix]=MC1.subs({x:xC1Num[ix]})
#
#for ix in range(Nx):
#    uC2Num[ix]=U_C.subs({x:xC2Num[ix]})
#    vC2Num[ix]=V_C2.subs({x:xC2Num[ix]})
#    PCNum[ix] = PC.subs({x:xC1Num[ix]})
#    VC2Num[ix]=VC2.subs({x:xC1Num[ix]})
#    MC2Num[ix]=MC2.subs({x:xC1Num[ix]})
    
for ix in range(Nx):
    if xCNum[ix] < float(LC/2):
        uCNum[ix]=U_C.subs({x:xCNum[ix]})
        vCNum[ix]=V_C1.subs({x:xCNum[ix]})
        PCNum[ix] = PC.subs({x:xCNum[ix]})
        VCNum[ix]=VC1.subs({x:xCNum[ix]})
        MCNum[ix]=MC1.subs({x:xCNum[ix]})
    
    else:
        uCNum[ix]=U_C.subs({x:xCNum[ix]})
        vCNum[ix]=V_C2.subs({x:xCNum[ix]})
        PCNum[ix] = PC.subs({x:xCNum[ix]})
        VCNum[ix]=VC2.subs({x:xCNum[ix]})
        MCNum[ix]=MC2.subs({x:xCNum[ix]})
    
    
        

for ix in range(Nx):
    uDNum[ix]=U_D.subs({x:xDNum[ix]})
    vDNum[ix]=V_D.subs({x:xDNum[ix]})
    PDNum[ix]=PD.subs({x:xDNum[ix]})
    VDNum[ix]=VD.subs({x:xDNum[ix]})
    MDNum[ix]=MD.subs({x:xDNum[ix]})

for ix in range(Nx):
    uENum[ix]=U_E.subs({x:xENum[ix]})
    vENum[ix]=V_E.subs({x:xENum[ix]})
    PENum[ix]=PE.subs({x:xENum[ix]})
    VENum[ix]=VE.subs({x:xENum[ix]})
    MENum[ix]=ME.subs({x:xENum[ix]})


xAFin=xANum+uANum*Escala
yAFin=0+vANum*Escala

xBFin=xBNum+uBNum*Escala
yBFin=0+vBNum*Escala

xCFin=xCNum+uCNum*Escala
yCFin=0+vCNum*Escala

#xC2Fin=xC2Num+uC2Num*Escala
#yC2Fin=0+vC2Num*Escala


xDFin=xDNum+uDNum*Escala
yDFin=0+vDNum*Escala

xEFin=xENum+uENum*Escala
yEFin=0+vENum*Escala

xAGloFig=np.cos(float(Theta_A))*xAFin-np.sin(float(Theta_A))*yAFin
yAGloFig=np.sin(float(Theta_A))*xAFin+np.cos(float(Theta_A))*yAFin

xBGloFig=np.cos(float(Theta_B))*xBFin-np.sin(float(Theta_B))*yBFin+x2
yBGloFig=np.sin(float(Theta_B))*xBFin+np.cos(float(Theta_B))*yBFin+y2

xCGloFig=np.cos(float(Theta_C))*xCFin-np.sin(float(Theta_C))*yCFin+x3
yCGloFig=np.sin(float(Theta_C))*xCFin+np.cos(float(Theta_C))*yCFin+y3

#xC2GloFig=np.cos(float(Theta_C))*xC2Fin-np.sin(float(Theta_C))*yC2Fin+x3
#yC2GloFig=np.sin(float(Theta_C))*xC2Fin+np.cos(float(Theta_C))*yC2Fin+y3

xDGloFig=np.cos(float(Theta_D))*xDFin-np.sin(float(Theta_D))*yDFin+x4
yDGloFig=np.sin(float(Theta_D))*xDFin+np.cos(float(Theta_D))*yDFin+y4

xEGloFig=np.cos(float(Theta_E))*xEFin-np.sin(float(Theta_E))*yEFin+x6
yEGloFig=np.sin(float(Theta_E))*xEFin+np.cos(float(Theta_E))*yEFin+y6



plt.figure(1,figsize=(8,5))

#   Figura elementos

plt.plot([x1,x2],[y1,y2],color='k')
plt.plot([x2,x3],[y2,y3],color='k')
plt.plot([x3,x4],[y3,y4],color='k')
plt.plot([x4,x5],[y4,y5],color='k')
plt.plot([x5,x6],[y5,y6],color='k')

#   Figura Deformada
plt.plot(xAGloFig,yAGloFig,color='r',linestyle='--')
plt.plot(xBGloFig,yBGloFig,color='r',linestyle='--')
plt.plot(xCGloFig,yCGloFig,color='r',linestyle='--')
#plt.plot(xC2GloFig,yC2GloFig,color='magenta',linestyle='--')
plt.plot(xDGloFig,yDGloFig,color='r',linestyle='--')
plt.plot(xEGloFig,yEGloFig,color='r',linestyle='--')


plt.grid('on')
plt.axis('equal')
plt.xlabel(r'$x[m]$',fontsize=16)
plt.ylabel(r'$y[m]$',fontsize=16)
plt.tick_params(labelsize=16)
plt.title('Configuración deformada')
#plt.text(x1+0.5  ,y1+0.5,np.round(float(MA.subs({x:0 })),3),fontsize=16)
#plt.text(x2-0.7  ,y2+0.5,np.round(float(MA.subs({x:LA})),3),fontsize=16)
#plt.text(x3-0.8  ,y3,np.round(float(MB.subs({x:LB})),3),fontsize=16)
#plt.text(x4-0.5  ,y4-1.5,np.round(float(MD.subs({x:0})),3),fontsize=16)
#plt.text(x5-0.5  ,y5,np.round(float(MD.subs({x:LD})),3),fontsize=16)
#plt.text(x5-0.5  ,y6+1.5,np.round(float(ME.subs({x:0})),3),fontsize=16)

plt.savefig('Figuras\ConfiguraciónDeformada.pdf')


plt.figure(2,figsize=(8,5))

#   Figura elementos

plt.plot([x1,x2],[y1,y2],color='k')
plt.plot([x2,x3],[y2,y3],color='k')
plt.plot([x3,x4],[y3,y4],color='k')
plt.plot([x4,x5],[y4,y5],color='k')
plt.plot([x5,x6],[y5,y6],color='k')

# Figura Axial

plt.plot(xANum*np.cos(float(Theta_A))-PANum*EscalaP*np.sin(float(Theta_A)),xANum*np.sin(float(Theta_A))+PANum*EscalaP*np.cos(float(Theta_A)),color='k',linestyle='--')
plt.plot(x2+xBNum*np.cos(float(Theta_B))-PBNum*EscalaP*np.sin(float(Theta_B)),y2+xBNum*np.sin(float(Theta_B))+PBNum*EscalaP*np.cos(float(Theta_B)),color='k',linestyle='--')

plt.plot(x3+xCNum*np.cos(float(Theta_C))-PCNum*EscalaP*np.sin(float(Theta_C)),y3+xCNum*np.sin(float(Theta_C))+PCNum*EscalaP*np.cos(float(Theta_C)),color='k',linestyle='--')
#plt.plot(x3+xC2Num*np.cos(float(Theta_C))-PCNum*EscalaP*np.sin(float(Theta_C)),y3+xC2Num*np.sin(float(Theta_C))+PCNum*EscalaP*np.cos(float(Theta_C)),color='k',linestyle='--')

plt.plot(x4+xDNum*np.cos(float(Theta_D))-PDNum*EscalaP*np.sin(float(Theta_D)),y4+xDNum*np.sin(float(Theta_D))+PDNum*EscalaP*np.cos(float(Theta_D)),color='k',linestyle='--')
plt.plot(x6+xENum*np.cos(float(Theta_E))-PENum*EscalaP*np.sin(float(Theta_E)),y6+xENum*np.sin(float(Theta_E))+PENum*EscalaP*np.cos(float(Theta_E)),color='k',linestyle='--')

plt.grid('on')
plt.axis('equal')
plt.xlabel(r'$x[m]$',fontsize=16)
plt.ylabel(r'$y[m]$',fontsize=16)
plt.tick_params(labelsize=16)
plt.title('Campo Fuerza Axial')
plt.text(x1+0.5  ,y1+0.5,np.round(float(PA.subs({x:0 })),3),fontsize=16)
plt.text(x2-0.7  ,y2+0.5,np.round(float(PA.subs({x:LA})),3),fontsize=16)
plt.text(x3-0.8  ,y3,np.round(float(PB.subs({x:LB})),3),fontsize=16)
plt.text(x4-0.5  ,y4-1.5,np.round(float(PC.subs({x:0})),3),fontsize=16)
plt.text(x5-2.5  ,y5,np.round(float(PD.subs({x:LD})),3),fontsize=16)
plt.text(x6-2.5  ,y6+0.3,np.round(float(PE.subs({x:0})),3),fontsize=16)

plt.savefig('Figuras\Campo Fuerza Axial.pdf')

plt.figure(3,figsize=(8,5))

#   Figura elementos

plt.plot([x1,x2],[y1,y2],color='k')
plt.plot([x2,x3],[y2,y3],color='k')
plt.plot([x3,x4],[y3,y4],color='k')
plt.plot([x4,x5],[y4,y5],color='k')
plt.plot([x5,x6],[y5,y6],color='k')

# Figura Cortante

plt.plot(xANum*np.cos(float(Theta_A))-VANum*EscalaV*np.sin(float(Theta_A)),xANum*np.sin(float(Theta_A))+VANum*EscalaV*np.cos(float(Theta_A)),color='k',linestyle='--')
plt.plot(x2+xBNum*np.cos(float(Theta_B))-VBNum*EscalaV*np.sin(float(Theta_B)),y2+xBNum*np.sin(float(Theta_B))+VBNum*EscalaV*np.cos(float(Theta_B)),color='k',linestyle='--')

plt.plot(x3+xCNum*np.cos(float(Theta_C))-VCNum*EscalaV*np.sin(float(Theta_C)),y3+xCNum*np.sin(float(Theta_C))+VCNum*EscalaV*np.cos(float(Theta_C)),color='k',linestyle='--')
#plt.plot(x3+xC2Num*np.cos(float(Theta_C))-VC2Num*EscalaV*np.sin(float(Theta_C)),y3+xC2Num*np.sin(float(Theta_C))+VC2Num*EscalaV*np.cos(float(Theta_C)),color='k',linestyle='--')

plt.plot(x4+xDNum*np.cos(float(Theta_D))-VDNum*EscalaV*np.sin(float(Theta_D)),y4+xDNum*np.sin(float(Theta_D))+VDNum*EscalaV*np.cos(float(Theta_D)),color='k',linestyle='--')
plt.plot(x6+xENum*np.cos(float(Theta_E))-VENum*EscalaV*np.sin(float(Theta_E)),y6+xENum*np.sin(float(Theta_E))+VENum*EscalaV*np.cos(float(Theta_E)),color='k',linestyle='--')

plt.grid('on')
plt.axis('equal')
plt.xlabel(r'$x[m]$',fontsize=16)
plt.ylabel(r'$y[m]$',fontsize=16)
plt.tick_params(labelsize=16)
plt.title('Campo Fuerza Cortante')
plt.text(x1+0.5  ,y1+0.5,np.round(float(VA.subs({x:0 })),3),fontsize=16)
plt.text(x2-0.7  ,y2+0.5,np.round(float(VA.subs({x:LA})),3),fontsize=16)
plt.text(x3-1.4  ,y3-0.3,np.round(float(VB.subs({x:LB})),3),fontsize=16)
plt.text(x4-1.2  ,y4-0.6,np.round(float(VD.subs({x:0})),3),fontsize=16)
plt.text(x5-2.5  ,y5,np.round(float(VD.subs({x:LD})),3),fontsize=16)
plt.text(x5-2.0  ,y6+0.3,np.round(float(VE.subs({x:0})),3),fontsize=16)

plt.savefig('Figuras\Campo Fuerza Cortante.pdf')

plt.figure(4,figsize=(8,5))

#   Figura elementos

plt.plot([x1,x2],[y1,y2],color='k')
plt.plot([x2,x3],[y2,y3],color='k')
plt.plot([x3,x4],[y3,y4],color='k')
plt.plot([x4,x5],[y4,y5],color='k')
plt.plot([x5,x6],[y5,y6],color='k')

# Figura Momentos

plt.plot(xANum*np.cos(float(Theta_A))-MANum*EscalaM*np.sin(float(Theta_A)),xANum*np.sin(float(Theta_A))+MANum*EscalaM*np.cos(float(Theta_A)),color='k',linestyle='--')
plt.plot(x2+xBNum*np.cos(float(Theta_B))-MBNum*EscalaM*np.sin(float(Theta_B)),y2+xBNum*np.sin(float(Theta_B))+MBNum*EscalaM*np.cos(float(Theta_B)),color='k',linestyle='--')

plt.plot(x3+xCNum*np.cos(float(Theta_C))-MCNum*EscalaM*np.sin(float(Theta_C)),y3+xCNum*np.sin(float(Theta_C))+MCNum*EscalaM*np.cos(float(Theta_C)),color='k',linestyle='--')
#plt.plot(x3+xC2Num*np.cos(float(Theta_C))-MC2Num*EscalaM*np.sin(float(Theta_C)),y3+xC2Num*np.sin(float(Theta_C))+MC2Num*EscalaM*np.cos(float(Theta_C)),color='k',linestyle='--')

plt.plot(x4+xDNum*np.cos(float(Theta_D))-MDNum*EscalaM*np.sin(float(Theta_D)),y4+xDNum*np.sin(float(Theta_D))+MDNum*EscalaM*np.cos(float(Theta_D)),color='k',linestyle='--')
plt.plot(x6+xENum*np.cos(float(Theta_E))-MENum*EscalaM*np.sin(float(Theta_E)),y6+xENum*np.sin(float(Theta_E))+MENum*EscalaM*np.cos(float(Theta_E)),color='k',linestyle='--')

plt.grid('on')
plt.axis('equal')
plt.xlabel(r'$x[m]$',fontsize=16)
plt.ylabel(r'$y[m]$',fontsize=16)
plt.tick_params(labelsize=16)
plt.title('Campo Momento Flector')
plt.text(x1+0.5  ,y1+0.5,np.round(float(MA.subs({x:0 })),3),fontsize=16)
plt.text(x2-0.7  ,y2+0.5,np.round(float(MA.subs({x:LA})),3),fontsize=16)
plt.text(x3-0.8  ,y3-0.8,np.round(float(MB.subs({x:LB})),3),fontsize=16)
plt.text(x4-0.5  ,y4-1.0,np.round(float(MD.subs({x:0})),3),fontsize=16)
plt.text(x5-2.0  ,y5,np.round(float(MD.subs({x:LD})),3),fontsize=16)
plt.text(x5-0.7  ,y6+0.3,np.round(float(ME.subs({x:0})),3),fontsize=16)

plt.savefig('Figuras\Campo Momento Flector.pdf')

#
##%% Funtions to latex
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(30)]
#    print(x)
#
#def funciones(F,n):
#    return sy.latex(sy.simplify(sy.N(F,n)))
#
#def Matrices(M,n):
#    ma = np.round(M,n)
#    sma = sy.Matrix(ma)
#    return sy.latex(sma)
#
#Nombres=['elemento A','elemento B','elemento C','elemento D','elemento E']
#
#Matrices_locales = [Matrices(Mloc_A,2),Matrices(Mloc_B,2),
#                    Matrices(Mloc_C,2),Matrices(Mloc_D,2),Matrices(Mloc_E,2)]
#
#Matrices_transformación = [Matrices(Mtrans_A,2),Matrices(Mtrans_B,2),
#                    Matrices(Mtrans_C,2),Matrices(Mtrans_D,2),Matrices(Mtrans_E,2)]
#
#fuerzas_empotramiento = [Matrices(Vemp_A,2),Matrices(Vemp_B,2),Matrices(Vemp_C,2),
#                         Matrices(Vemp_D,2),Matrices(Vemp_E,2)]
#
#Matrices_Globales = [Matrices(KA,2),Matrices(KB,2),Matrices(KC,2),Matrices(KD,2),Matrices(KE,2)]
#
#fuerzas_empotramiento_Glo = [Matrices(Vemp_GA,2),Matrices(Vemp_GB,2),Matrices(Vemp_GC,2),
#                         Matrices(Vemp_GD,2),Matrices(Vemp_GE,2)]
#
#Desplazamientos_locales=[Matrices(DesNod_Glo_A,2),Matrices(DesNod_Glo_B,2),Matrices(DesNod_Glo_C,2),
#                        Matrices(DesNod_Glo_D,2),Matrices(DesNod_Glo_E,2)]
#print('Matrices Locales \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(Matrices_locales[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('Matrices De transformación \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(Matrices_transformación[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('fuerzas de empotramiento en locales\n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(fuerzas_empotramiento[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('Matrices Globales \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(Matrices_Globales[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('fuerzas de empotramiento_Globales \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(fuerzas_empotramiento_Glo[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('Desplazamientos en coordenadas locales \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(Desplazamientos_locales[k])
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('MATRIZ GRANDE DE RIGIDEZ \n', Matrices(DesNod,2) ,'\n','Vector de Fuerzas de Empotramiento \n',Matrices(FUE_EMP,2))
#
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('Resultados con notación   \n',sy.latex(des_notacion))
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#
#print('Transformación desplzamientos a coordenadas locales \n')
#
#for k in range(len(Nombres)):
#    print(Nombres[k])
#    print(Desplazamientos_locales[k])
#    
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#    
#print('Calculo de las reacciones \n','FX1=',np.round(FX1,3),'\n','FY1=',np.round(FY1,3),'\n',
#      'FX6=',np.round(FX6,3),'\n','FY6=',np.round(FY6,3),'\n')
#
#KREA = np.zeros([4,15])
#
#KREA[0] = DesNod[1]
#KREA[1] = DesNod[2]
#KREA[2] = DesNod[11]
#KREA[3] = DesNod[12]
#
#print('Matriz porción reacciones \n', Matrices(KREA,2),'\n')
#
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)
#    
#nombres_campos = ['U_E','V_E']
#nombres_campos_C = ['U_C','V_C1','V_C2']
#Campos_des_loc_A = [funciones(U_A,3),funciones(V_A,3)]
#Campos_des_loc_B = [funciones(U_B,3),funciones(V_B,3)]
#Campos_des_loc_C = [funciones(U_C,3),funciones(V_C1,3),funciones(V_C2,3)]
#Campos_des_loc_D = [funciones(U_D,3),funciones(V_D,3)]
#Campos_des_loc_E = [funciones(U_E,3),funciones(V_E,3)]
#
#Fuerza_suelo = [funciones(F_soil_A,3),funciones(F_soil_E,3)]
#
#fuerzas = ['PE','VA','MA']
#fuerzas_C = ['PC','VC1','VC2','MC1','MC2']
#
#Fuerzas_internas_A = [funciones(PA,3),funciones(VA,3),funciones(MA,3)]
#Fuerzas_internas_B = [funciones(PB,3),funciones(VB,3),funciones(MB,3)]
#Fuerzas_internas_C = [funciones(PC,3),funciones(VC1,3),funciones(VC2,3),funciones(MC1,3),funciones(MC2,3)]
#Fuerzas_internas_D = [funciones(PD,3),funciones(VD,3),funciones(MD,3)]
#Fuerzas_internas_E = [funciones(PE,3),funciones(VE,3),funciones(ME,3)]
#
#print('Campos de Desplazameintos en locales  \n')
#print('Elemento A \n')
#for kl in range(len(nombres_campos)):
#    print(nombres_campos[kl],'\n',(Campos_des_loc_A)[kl],'\n', '\n')
#print('Elemento B \n')    
#for kl in range(len(nombres_campos)):
#    print(nombres_campos[kl],'\n',(Campos_des_loc_B)[kl],'\n', '\n')
#print('Elemento C \n')
#for kl in range(len(nombres_campos_C)):
#    print(nombres_campos_C[kl],'\n',(Campos_des_loc_C)[kl],'\n', '\n')
#print('Elemento D \n')
#for kl in range(len(nombres_campos)):
#    print(nombres_campos[kl],'\n',(Campos_des_loc_D)[kl],'\n', '\n')
#print('Elemento E \n')
#for kl in range(len(nombres_campos)):
#    print(nombres_campos[kl],'\n',(Campos_des_loc_E)[kl],'\n', '\n')
#        
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)    
#
#print('Fuerza de suelos sobre las pilas \n ')
#
##for i in range(len(Fuerza_suelo)):
##    if i == 0:
##        print('\underline{Elemento A} \n',Fuerza_suelo[i],'\n')
##    else:
##        print('\underline{Elemento E} \n',Fuerza_suelo[i],'\n')
#    
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)    
#print('\n')
#for ll in range(len(fuerzas)):
#    print(Nombres[0],'\n')
#    print(fuerzas[ll],'\n',Fuerzas_internas_A[ll],'\n', '\n')
#    
#for ll in range(len(fuerzas)):
#    print(Nombres[1],'\n')
#    print(fuerzas[ll],'\n',Fuerzas_internas_B[ll],'\n', '\n')
#
#for ll in range(len(fuerzas)):
#    print(Nombres[2],'\n')
#    print(fuerzas_C[ll],'\n',Fuerzas_internas_C[ll],'\n', '\n')
#
#for ll in range(len(fuerzas)):
#    print(Nombres[3],'\n')
#    print(fuerzas[ll],'\n',Fuerzas_internas_D[ll],'\n', '\n')
#    
#for ll in range(len(fuerzas)):
#    print(Nombres[4],'\n')
#    print(fuerzas[ll],'\n',Fuerzas_internas_E[ll],'\n', '\n')
#    
#for ij in range(1):
#    ls=[]
#    x = ['#' for i in range(20)]
#    print(x)   

