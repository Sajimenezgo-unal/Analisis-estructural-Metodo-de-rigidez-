# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 16:54:12 2021

@author: Sajimenezgo 
"""
### Trabajo Estructural ###

# Librerias

import sympy as sy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import init_printing
init_printing() 


#%% funciones 

def radians(x):
    r = (x*sy.pi)/180
    return r

##          Funciones de funciones de forma 

def Fforma(x,L):
    Psi1 = 1 - x/L
    Psi2 = 1 - 3*(x/L)**2 + 2*(x/L)**3
    Psi3 = ((x/L)-2*(x/L)**2 + (x/L)**3)*L
    Psi4 = x/L
    Psi5 = 3*(x/L)**2 -2*(x/L)**3
    Psi6 = (-1*(x/L)**2 + (x/L)**3)*L
    
    list = [Psi1,Psi2,Psi3,Psi4,Psi5,Psi6]
    return list

##          Funciones de green 

def Gyy(x,xi,Le,Ee,Ie):
    gyy1 = ((Le**3)/(6*Ee*Ie))*(-1*(x/Le)**3*Fforma(xi,Le)[1]
                                + (3*(x/Le)**2)*Fforma(xi,Le)[2]/Le)
    gyy2 = ((Le**3)/(6*Ee*Ie))*(-1*(1- (x/Le))**3*Fforma(xi,Le)[4]
                                - (3*(1- (x/Le))**2)*Fforma(xi,Le)[5]/Le)
    list = [gyy1,gyy2]
    return list 

def Gyr(x,xi,Le,Ee,Ie):
    gyr1 = ((Le**2)/(6*Ee*Ie))*((-1*(x/Le)**3)*(sy.diff(Fforma(xi,Le)[1],xi))
                                + (3*(x/Le)**2)*(sy.diff(Fforma(xi,Le)[2],xi))/Le)
    gyr2 = ((Le**2)/(6*Ee*Ie))*(-1*(1- (x/Le))**3*(sy.diff(Fforma(xi,Le)[4],xi))
                                - (3*(1- (x/Le))**2)*(sy.diff(Fforma(xi,Le)[5],xi))/Le)
    list = [gyr1,gyr2]
    return list 

def Gxx(x,Le,Ee,Ae):
    gxx1 = (Le/Ae*Ee)*Fforma(x,Le)[3]*Fforma(x,Le)[0]
    gxx2 = (Le/Ae*Ee)*Fforma(x,Le)[0]*Fforma(x,Le)[3]
    list = [gxx1,gxx2]
    return list 

### Elemento tipo portico plano 

def cdeshomogeneo(ui,uj,vi,vj,Theta_i,Theta_j,x,L):
    Ueh = Fforma(x,L)[0]*ui + Fforma(x,L)[3]*uj
    Veh = Fforma(x,L)[1]*vi + Fforma(x,L)[2]*Theta_i + Fforma(x,L)[4]*vj + Fforma(x,L)[5]*Theta_j
    l = [Ueh,Veh]
    return l

def cdesempotrado(xi,xe,x,Le,Ee,Ae,p,q):
    
    Uef = sy.integrate(p.subs({x:xi})*Gxx(x,Le,Ee,Ae)[1],(xi,0,xe)) + sy.integrate(p.subs({x:xi})*Gxx(x,Le,Ee,Ae)[0],(xi,xe,Le))
    
    Vef = sy.integrate(q.subs({x:xi})*Gyy(x,xi,Le,Ee,Ie)[1],(xi,0,xe)) + sy.integrate(q.subs({x:xi})*Gyy(x,xi,Le,Ee,Ie)[0],(xi,xe,Le))
    
    l = [Uef,Vef]
    return l

def Matrizportico_locales(Ae,Ee,Le,Ie):
    
    MLoc = sy.zeros(6,6)
    
    MLoc[0,0]= (Ae * Ee)/Le
    MLoc[0,1]=0
    MLoc[0,2]=0
    MLoc[0,3]= -1*(Ae * Ee)/Le              #fila 1
    MLoc[0,4]=0
    MLoc[0,5]=0
    
    MLoc[1,0]=0
    MLoc[1,1]= 12*Ee*Ie/(Le**3)
    MLoc[1,2]= 6*Ee*Ie/(Le**2)
    MLoc[1,3]= 0                            #fila 2
    MLoc[1,4]= -12*Ee*Ie/(Le**3)
    MLoc[1,5]= 6*Ee*Ie/(Le**2)

    MLoc[2,0]=0
    MLoc[2,1]= 6*Ee*Ie/(Le**2)
    MLoc[2,2]= 4*Ee*Ie/(Le)
    MLoc[2,3]= 0                            #fila 3
    MLoc[2,4]= -6*Ee*Ie/(Le**2)
    MLoc[2,5]= 2*Ee*Ie/(Le)

    MLoc[3,0]=-1*(Ae * Ee)/Le
    MLoc[3,1]=0
    MLoc[3,2]=0
    MLoc[3,3]=(Ae * Ee)/Le                  #fila 4
    MLoc[3,4]=0
    MLoc[3,5]=0

    MLoc[4,0]=0
    MLoc[4,1]= -12*Ee*Ie/(Le**3)
    MLoc[4,2]= -6*Ee*Ie/(Le**2)
    MLoc[4,3]=0                             #fila 5
    MLoc[4,4]= 12*Ee*Ie/(Le**3)
    MLoc[4,5]= -6*Ee*Ie/(Le**2)

    MLoc[5,0]= 0
    MLoc[5,1]= 6*Ee*Ie/(Le**2)
    MLoc[5,2]= 2*Ee*Ie/(Le)
    MLoc[5,3]= 0                            #fila 6
    MLoc[5,4]= -6*Ee*Ie/(Le**2)
    MLoc[5,5]= 4*Ee*Ie/(Le)
    
    return MLoc

def Vector_emp_portico(x,L,p,q):
    vec = sy.zeros(6,1, dtype = float)
    
    vec[0,0] = -sy.integrate(Fforma(x,L)[0]*p,(x,0,L))
    vec[1,0] = -sy.integrate(Fforma(x,L)[1]*q,(x,0,L))
    vec[2,0] = -sy.integrate(Fforma(x,L)[2]*q,(x,0,L))
    vec[3,0] = -sy.integrate(Fforma(x,L)[3]*p,(x,0,L))
    vec[4,0] = -sy.integrate(Fforma(x,L)[4]*q,(x,0,L))
    vec[5,0] = -sy.integrate(Fforma(x,L)[5]*q,(x,0,L))
    

    
    return vec

def Mtrans(theta_e,phi_i,phi_j):
    TE = sy.zeros(6,6)
    
    TE[0,0] = sy.cos(theta_e - phi_i)
    TE[0,1] = sy.sin(theta_e - phi_i)
    
    TE[1,0] = -1*sy.sin(theta_e - phi_i)
    TE[1,1] = sy.cos(theta_e-phi_i)

    TE[2,2] = 1

    TE[3,3] = sy.cos(theta_e-phi_j)
    TE[3,4] = sy.sin(theta_e-phi_j)

    TE[4,3] = -1*sy.sin(theta_e-phi_j)
    TE[4,4] = sy.cos(theta_e-phi_j)

    TE[5,5] = 1 
    return TE

## Elemento tipo viga sobre fundación flexible 

def Fforma_pila(ke,Le,Ie,Ee,xe):

    lambda_e = (ke/(4*Ee*Ie))**(1/4)
    s = sy.sin(lambda_e*Le)
    c = sy.cos(lambda_e*Le)
    sh = sy.sinh(lambda_e*Le)
    ch = sy.cosh(lambda_e*Le)
    
    psi2 = (-1*((s*ch)**2 + (c*sh)**2)*sy.sin(lambda_e*xe)*sy.sinh(lambda_e*xe) + (s*c+sh*ch)*(sy.sin(lambda_e*xe)*sy.cosh(lambda_e*xe) - sy.cos(lambda_e*xe)*sy.sinh(lambda_e*xe)) +(sh**2-s**2)*sy.cos(lambda_e*xe)*sy.cosh(lambda_e*xe))/(sh**2-s**2)
    psi3 = psi3 = (1/lambda_e)*(((s*c-sh*ch)*sy.sin(lambda_e*xe)*sy.sinh(lambda_e*xe)+sh**2*sy.sin(lambda_e*xe)*sy.cosh(lambda_e*xe)-s**2*sy.cos(lambda_e*xe)*sy.sinh(lambda_e*xe))/(sh**2 - s**2))
    psi5 = (2*s*sh*sy.sin(lambda_e*xe)*sy.sinh(lambda_e*xe)-(s*ch+c*sh)*sy.sin(lambda_e*xe)*sy.cosh(lambda_e*xe)+(s*ch+c*sh)*sy.cos(lambda_e*xe)*sy.sinh(lambda_e*xe))/(sh**2 -s**2)
    psi6 = (1/lambda_e)*(((c*sh-s*ch)*sy.sin(lambda_e*xe)*sy.sinh(lambda_e*xe)+s*sh*sy.sin(lambda_e*xe)*sy.cosh(lambda_e*xe)-s*sh*sy.cos(lambda_e*xe)*sy.sinh(lambda_e*xe))/(sh**2 - s**2))
    l = [0,psi2,psi3,0,psi5,psi6]
    return l

def gyy_pilas(ke,Le,Ie,Ee,xe,xi):
    lambda_e = (ke/(4*Ee*Ie))**(1/4)
    gyy1 = (1/Ee*Ie)*(((-1*(sy.sin(lambda_e*xe)*sy.cosh(lambda_e*xe))-sy.sinh(lambda_e*xe)*sy.cos(lambda_e*xe))/(4*lambda_e**3))*(Fforma_pila(ke,Le,Ie,Ee,xi)[1]) + ((sy.sin(lambda_e*xe)*sy.sinh(lambda_e*xe))/(2*lambda_e**2))*(Fforma_pila(ke,Le,Ie,Ee,xi)[2]))
    gyy2 = (1/Ee*Ie)*(((-1*(sy.sin(lambda_e*(Le-xe))*sy.cosh(lambda_e*(Le-xe)))-sy.sinh(lambda_e*(Le-xe))*sy.cos(lambda_e*(Le-xe)))/(4*lambda_e**3))*(Fforma_pila(ke,Le,Ie,Ee,xi)[4]) - ((sy.sin(lambda_e*(Le-xe))*sy.sinh(lambda_e*(Le-xe)))/(2*lambda_e**2))*(Fforma_pila(ke,Le,Ie,Ee,xi)[5]))
    
    l = [gyy1,gyy2]
    return l


def Matriz_pilas_locales(ke,Le,Ie,Ee,xe,Ae):
    
    lambda_e = (ke/(4*Ee*Ie))**(1/4)
    s = sy.sin(lambda_e*Le)
    c = sy.cos(lambda_e*Le)
    sh = sy.sinh(lambda_e*Le)
    ch = sy.cosh(lambda_e*Le)
    
    M = sy.zeros(6,6, dtype = float)
    k22=k55 = 4*Ee*Ie*(lambda_e**3)*((s*c + sh*ch)/((sh**2) - s**2))
    k23=k32 = 2*Ee*Ie*(lambda_e**2)*(((s**2) + sh**2)/((sh**2) - s**2))
    k56=k65 = -2*Ee*Ie*(lambda_e**2)*(((s**2) + sh**2)/((sh**2) - s**2))
    k25=k52 = -4*Ee*Ie*(lambda_e**3)*((s*ch + c*sh)/((sh**2)-s**2))
    k26=k62 = 4*Ee*Ie*(lambda_e**2)*((s*sh)/((sh**2)-s**2))
    k35=k53 = -4*Ee*Ie*(lambda_e**2)*((s*sh)/((sh**2)-s**2))
    k33=k66 = 2*Ee*Ie*(lambda_e)*((sh*ch - s*c)/((sh**2) - s**2))
    k36=k63 = 2*Ee*Ie*(lambda_e)*((s*ch - c*sh)/((sh**2) - s**2))
    
    # fila 1
    M[0,0] = Ae*Ee/Le 
    M[0,1] = 0
    M[0,2] = 0
    M[0,3] = -Ae*Ee/Le
    M[0,4] = 0
    M[0,5] = 0 

    # fila 2
    M[1,0] = 0 
    M[1,1] = k22
    M[1,2] = k23
    M[1,3] = 0 
    M[1,4] = k25
    M[1,5] = k26
    
    # fila 3
    M[2,0] = 0 
    M[2,1] = k32
    M[2,2] = k33
    M[2,3] = 0
    M[2,4] = k35
    M[2,5] = k36
    
    # fila 4
    M[3,0] = -Ae*Ee/Le
    M[3,1] = 0
    M[3,2] = 0 
    M[3,3] = Ae*Ee/Le
    M[3,4] = 0
    M[3,5] = 0

    # fila 5
    M[4,0] = 0
    M[4,1] = k52
    M[4,2] = k53
    M[4,3] = 0
    M[4,4] = k55
    M[4,5] = k56
    
    # fila 6
    M[5,0] = 0 
    M[5,1] = k62
    M[5,2] = k63
    M[5,3] = 0
    M[5,4] = k65
    M[5,5] = k66
    
    return M

def vector_em_pilas(ke,Le,Ie,Ee,xe,q,p):
    vect = sy.zeros(6,1, dtype = float)
    
    vect[0,0] = -sy.integrate(Fforma(xe,Le)[0]*p,(xe,0,Le))
    vect[1,0] = -sy.integrate(Fforma_pila(ke,Le,Ie,Ee,xe)[1]*q,(xe,0,Le))
    vect[2,0] = -sy.integrate(Fforma_pila(ke,Le,Ie,Ee,xe)[2]*q,(xe,0,Le))
    vect[3,0] = -sy.integrate(Fforma(xe,Le)[3]*p,(xe,0,Le))
    vect[4,0] = -sy.integrate(Fforma_pila(ke,Le,Ie,Ee,xe)[4]*q,(xe,0,Le))
    vect[5,0] = -sy.integrate(Fforma_pila(ke,Le,Ie,Ee,xe)[5]*q,(xe,0,Le))
    
    return vect
    
def Cargas2Locales (dirR,Q,theta):
    if dirR == 'v' or 'V':
        p = Q*sy.Abs(sy.cos(theta))*sy.sin(theta)
        q = Q*sy.Abs(sy.cos(theta))*sy.cos(theta)
        
    else:
        p = Q*sy.Abs(sy.sin(theta))*sy.cos(theta)
        q = -1*Q*sy.Abs(sy.sin(theta))*sy.sin(theta)
        
    l = [p,q]
    return l
#%% Matrices 
xe,xi,Le,Ee,Ie,p,q = sy.symbols('xe,xi,Le,Ee,Ie,p,q')
ui,uj,vi,vj,thetai,thetaj = sy.symbols('ui,uj,vi,vj,Theta_i,Theta_j')
Ae,Ee,Le,Ie,phi_i,phi_j,theta_e = sy.symbols('Ae,Ee,Le,Ie,phi_i,phi_j,theta_e')
# A = Fforma(xa,La)[0]
Le,Ie,Ee,xe = sy.symbols('Le,Ie,Ee,xe')
ke,Le = sy.symbols('ke,Le', positive = True)

#%% Información Elementos

LA = 5                          #longitud de A en metros 
LB = 4                          #longitud de B en metros
LC = sy.sqrt(6**2 + 2**2)       #longitud de C en metros
LD = sy.sqrt(5**2 + 2**2)       #longitud de D en metros

E = 2.0*(10**7)

r_pila = (110/100)/2
Area_pila = sy.pi*(r_pila**2)

b,h = 0.40,0.30
Area_portico = b*h

k = 9600
 
#%% transformación de cargas de globales a locales

x,xa,xb,xc,xd = sy.symbols('x,xa,xb,xc,xd')

############### Cargas Distribuidas en los Elementos Globales

Qc = (10/3)*xc - 90 ##  3m< x <6m      ## Hacia arriba 
Qd = (70/5)*xd - 70 ##  0m< x <5m     ## Hacia arriba

############### Cargas Distribuidas en los Elementos Locales

theta_C = sy.atan2((2),(6))
theta_D = sy.atan2((-2),(5))

p_C = Cargas2Locales ('V',Qc,theta_C)[0]        ### Elemento C
q_C = Cargas2Locales ('V',Qc,theta_C)[1]        ### Elemento C

p_D = Cargas2Locales ('V',Qd,theta_D)[0]        ### Elemento D
q_D = Cargas2Locales ('V',Qd,theta_D)[1]        ### Elemento D

#%% Definición de los modelos matriciales en coordenadas locales

### Elemento tipo pila ---- A

Ix_pila = (1/4)*sy.pi*(r_pila**4) #Inercia para elemento tipo pila
MA_loc = Matriz_pilas_locales(k,LA,Ix_pila,E,xa,Area_pila)
Vec_Emp_A =  vector_em_pilas(k,LA,Ix_pila,E,xa,0,0)



### Elemento tipo portico ---- B,C,D

Ix_P = (1/12)*b*(h**3) #Inercia para elementos tipo portico

MB_loc = Matrizportico_locales(Area_portico,E,LB,Ix_P)
Vec_Emp_B = Vector_emp_portico(xb,LB,0,0)

MC_loc = Matrizportico_locales(Area_portico,E,LC,Ix_P)
Vec_Emp_C = Vector_emp_portico(xc,LC,p_C,sy.simplify(q_C, rational =True))

MD_loc = Matrizportico_locales(Area_portico,E,LD,Ix_P)
Vec_Emp_D = Vector_emp_portico(xd,LD,p_D,sy.simplify(q_D, rational =True))


#%% Transformación de matrices a coordenadas globales

MAtrans = Mtrans(radians(90),0,0)
MBtrans = Mtrans(radians(90),0,0)
MCtrans = Mtrans(theta_C,0,0)
MDtrans = Mtrans(theta_D,0,0)

KA = (np.transpose(MAtrans))@ MA_loc @ MAtrans
KB = (np.transpose(MBtrans))@ MB_loc @ MBtrans
KC = (np.transpose(MCtrans))@ MC_loc @ MCtrans
KD = (np.transpose(MDtrans))@ MD_loc @ MDtrans

VA_Emp_Glo = (np.transpose(MAtrans))@ Vec_Emp_A 
VB_Emp_Glo = (np.transpose(MBtrans))@ Vec_Emp_B 
VC_Emp_Glo = (np.transpose(MCtrans))@ Vec_Emp_C 
VD_Emp_Glo = (np.transpose(MDtrans))@ Vec_Emp_D 

#%% Sistema Matricial 
Desnod = sy.zeros(7,7)

for i in range(3):
    for j in range(3):
        '''
        Suma de los componente de KA y KB de la siguiente manera 
        KA[:,3:] y KB[:,:3] para los desplazamientos de U2, V2, theta2 
        
        Suma de los componente de KB, KC y KD de la siguiente manera 
        KB[3:,3:], KC[3:,3:] y KD[:,:3] para los desplazamientos de U3, V3, theta3 
        
        '''
        Desnod[i,j] = KA[i,j+3] + KB[i,j]
        Desnod[i+3,j+3] = KB[i+3,j+3] + KC[i+3,j+3] + KD[i,j]
        
        # Verificación de la suma # 
        """ Comentar las lineas cuando no se necesiten """ 
        
        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KA[{0},{1}+3]'.format(i, j) + ' + ' + 'KB[{0},{1}]'.format(i, j))
        # print('\n' for i in range(5))
        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KB[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0},{1}]'.format(i, j))
        
Desnod[6,6] = KC[2,2]
# print(Desnod)

Vect_emp = sy.zeros(7,1)

for i in range(3):
    '''
    Suma de las componentes de los vectores de empotramiento 
    '''
    Vect_emp[i,0] = VA_Emp_Glo[i+3,0] + VB_Emp_Glo[i,0] 
    Vect_emp[i+3,0] = VB_Emp_Glo[i+3,0] + VC_Emp_Glo[i+3,0] + VD_Emp_Glo[i,0]
    
        # Verificación de la suma # 
    """ Comentar las lineas cuando no se necesiten """ 
        
    # print('VA_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VB_Emp_Glo[{},0]'.format(i)) 
    # print('VB_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VC_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VC_Emp_Glo[{},0]'.format(i))
    
#%% Solución del Sistema matricial 

"""  
Conversión a un array de pandas para utilizar la función linalg.solve()
"""

sis_Eq =  np.array(Desnod).astype(np.float64)
vec_emp =  np.array(-1*Vect_emp).astype(np.float64) 

Desplaza = np.linalg.solve(sis_Eq, vec_emp )

######## Desplazamientos en notacion cientifica ####################
l = []
for i in range(7):
    v = np.format_float_scientific(Desplaza[i,0], precision=2,exp_digits=1)
    l.append(v)
    
des_notacion = (sy.Matrix(np.transpose(np.asmatrix(l))))

####################################################################



#%% ensayos


Vect = sy.zeros(7,1)


for i in range(3):
	Vect[i,0] = VA_Emp_Glo[i+3,0] + VC_Emp_Glo[i,0]
	Vect[i+3,0] = VB_Emp_Glo[i+3,0] + VC_Emp_Glo[i+3,0] + VD_Emp_Glo[i,0]
Vect[6,0] = VC_Emp_Glo[2,0]
    
  

            