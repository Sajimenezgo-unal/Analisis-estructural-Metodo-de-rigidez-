# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 16:54:12 2021

@author: Sajimenezgo 
"""
### Trabajo Estructural ###

# Librerias


from Modulos import *
import sympy as sy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import init_printing
init_printing()

"""
IMPORTACIÓN DE LOS MODULOS LOS CUALES CONTIENEN TODAS LAS FORMULAS PARA EL DESARROLLO DEL EJERCICIO 

"""

# %% Matrices
xe, xi, Le, Ee, Ie, p, q = sy.symbols('xe,xi,Le,Ee,Ie,p,q')
ui, uj, vi, vj, thetai, thetaj = sy.symbols('ui,uj,vi,vj,Theta_i,Theta_j')
Ae, Ee, Le, Ie, phi_i, phi_j, theta_e = sy.symbols(
    'Ae,Ee,Le,Ie,phi_i,phi_j,theta_e')
# A = Fforma(xa,La)[0]
Le, Ie, Ee, xe = sy.symbols('Le,Ie,Ee,xe')
ke, Le = sy.symbols('ke,Le', positive=True)

# %% Información Elementos

LA = 5  # longitud de A en metros
LB = 4  # longitud de B en metros
LC = sy.sqrt(6**2 + 2**2)  # longitud de C en metros
LD = sy.sqrt(5**2 + 2**2)  # longitud de D en metros

E = 2.0*(10**7)

r_pila = (110/100)/2
Area_pila = sy.pi*(r_pila**2)

b, h = 0.40, 0.30
Area_portico = b*h

k = 9600

# %% transformación de cargas de globales a locales

x, xa, xb, xc, xd = sy.symbols('x,xa,xb,xc,xd')

# Cargas Distribuidas en los Elementos Globales

Qc = (10/3)*xc - 90  # 3m< x <6m      ## Hacia arriba
Qd = (70/5)*xd - 70  # 0m< x <5m     ## Hacia arriba

# Cargas Distribuidas en los Elementos Locales

theta_C = sy.atan2((2), (6))
theta_D = sy.atan2((-2), (5))

p_C = Cargas2Locales('V', Qc, theta_C)[0]  # Elemento C
q_C = Cargas2Locales('V', Qc, theta_C)[1]  # Elemento C

p_D = Cargas2Locales('V', Qd, theta_D)[0]  # Elemento D
q_D = Cargas2Locales('V', Qd, theta_D)[1]  # Elemento D

# %% Definición de los modelos matriciales en coordenadas locales

# Elemento tipo pila ---- A

Ix_pila = (1/4)*sy.pi*(r_pila**4)  # Inercia para elemento tipo pila
MA_loc = Matriz_pilas_locales(k, LA, Ix_pila, E, xa, Area_pila)
Vec_Emp_A = vector_em_pilas(k, LA, Ix_pila, E, xa, 0, 0)


# Elemento tipo portico ---- B,C,D

Ix_P = (1/12)*b*(h**3)  # Inercia para elementos tipo portico

MB_loc = Matrizportico_locales(Area_portico, E, LB, Ix_P)
Vec_Emp_B = Vector_emp_portico(xb, LB, 0, 0)

MC_loc = Matrizportico_locales(Area_portico, E, LC, Ix_P)
Vec_Emp_C = Vector_emp_portico(xc, LC, p_C, sy.simplify(q_C, rational=True))

MD_loc = Matrizportico_locales(Area_portico, E, LD, Ix_P)
Vec_Emp_D = Vector_emp_portico(xd, LD, p_D, sy.simplify(q_D, rational=True))


# %% Transformación de matrices a coordenadas globales

MAtrans = Mtrans(radians(90), 0, 0)
MBtrans = Mtrans(radians(90), 0, 0)
MCtrans = Mtrans(theta_C, 0, 0)
MDtrans = Mtrans(theta_D, 0, 0)

KA = (np.transpose(MAtrans)) @ MA_loc @ MAtrans
KB = (np.transpose(MBtrans)) @ MB_loc @ MBtrans
KC = (np.transpose(MCtrans)) @ MC_loc @ MCtrans
KD = (np.transpose(MDtrans)) @ MD_loc @ MDtrans

VA_Emp_Glo = (np.transpose(MAtrans)) @ Vec_Emp_A
VB_Emp_Glo = (np.transpose(MBtrans)) @ Vec_Emp_B
VC_Emp_Glo = (np.transpose(MCtrans)) @ Vec_Emp_C
VD_Emp_Glo = (np.transpose(MDtrans)) @ Vec_Emp_D

# %% Sistema Matricial
Desnod = sy.zeros(7, 7)

for i in range(3):
    for j in range(3):
        '''
        Suma de los componente de KA y KB de la siguiente manera 
        KA[:,3:] y KB[:,:3] para los desplazamientos de U2, V2, theta2 

        Suma de los componente de KB, KC y KD de la siguiente manera 
        KB[3:,3:], KC[3:,3:] y KD[:,:3] para los desplazamientos de U3, V3, theta3 

        '''
        Desnod[i, j] = KA[i, j+3] + KB[i, j]
        Desnod[i+3, j+3] = KB[i+3, j+3] + KC[i+3, j+3] + KD[i, j]

        # Verificación de la suma #
        """ Comentar las lineas cuando no se necesiten """

        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KA[{0},{1}+3]'.format(i, j) + ' + ' + 'KB[{0},{1}]'.format(i, j))
        # print('\n' for i in range(5))
        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KB[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0},{1}]'.format(i, j))

Desnod[6, 6] = KC[2, 2]
# print(Desnod)

Vect_emp = sy.zeros(7, 1)

for i in range(3):
    '''
    Suma de las componentes de los vectores de empotramiento 
    '''
    Vect_emp[i, 0] = VA_Emp_Glo[i+3, 0] + VB_Emp_Glo[i, 0]
    Vect_emp[i+3, 0] = VB_Emp_Glo[i+3, 0] + \
        VC_Emp_Glo[i+3, 0] + VD_Emp_Glo[i, 0]

    # Verificación de la suma #
    """ Comentar las lineas cuando no se necesiten """

    # print('VA_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VB_Emp_Glo[{},0]'.format(i))
    # print('VB_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VC_Emp_Glo[{}+3,0]'.format(i) + ' + ' + 'VC_Emp_Glo[{},0]'.format(i))

# %% Solución del Sistema matricial

"""  
Conversión a un array de pandas para utilizar la función linalg.solve()
"""

sis_Eq = np.array(Desnod).astype(np.float64)
vec_emp = np.array(-1*Vect_emp).astype(np.float64)

Desplaza = np.linalg.solve(sis_Eq, vec_emp)

######## Desplazamientos en notacion cientifica ####################
l = []
for i in range(7):
    v = np.format_float_scientific(Desplaza[i, 0], precision=2, exp_digits=1)
    l.append(v)

des_notacion = (sy.Matrix(np.transpose(np.asmatrix(l))))

####################################################################


# %% ensayos


Vect = sy.zeros(7, 1)


for i in range(3):
    Vect[i, 0] = VA_Emp_Glo[i+3, 0] + VC_Emp_Glo[i, 0]
    Vect[i+3, 0] = VB_Emp_Glo[i+3, 0] + VC_Emp_Glo[i+3, 0] + VD_Emp_Glo[i, 0]
Vect[6, 0] = VC_Emp_Glo[2, 0]
