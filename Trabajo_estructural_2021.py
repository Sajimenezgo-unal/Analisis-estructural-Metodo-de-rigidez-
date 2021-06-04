# %%

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
        Desnod[i, j] = KA[i+3, j+3] + KB[i, j]
        Desnod[i+3, j+3] = KB[i+3, j+3] + KC[i+3, j+3] + KD[i, j]
        Desnod[i, j+3] = KB[i, j+3]
        Desnod[i+3, j] = KB[i+3, j]
        Desnod[i+3, 6] = KC[i+3, 2]
        Desnod[6, j+3] = KC[2, j+3]
        # Verificación de la suma #
        """ Comentar las lineas cuando no se necesiten """

        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KA[{0},{1}+3]'.format(i, j) + ' + ' + 'KB[{0},{1}]'.format(i, j))
        # print('\n' for i in range(5))
        # print('Verificación de suma de los desplazamientos en el Nodo 2 \n ')
        # print('KB[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0}+3,{1}+3]'.format(i, j) + ' + ' + 'KC[{0},{1}]'.format(i, j))

Desnod[6, 6] = KC[2, 2]
# Desnod[6, 5] = KC[2, 5]
# Desnod[6, 4] = KC[2, 4]
# Desnod[6, 3] = KC[2, 3]
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

'''
Desplazamientos nodales en coordenadas globales
'''

U1 = 0
V1 = 0
Theta1 = 0

U2 = Desplaza[0, 0]
V2 = Desplaza[1, 0]
Theta2 = Desplaza[2, 0]

U3 = Desplaza[3, 0]
V3 = Desplaza[4, 0]
Theta3 = Desplaza[5, 0]


U4 = 0
V4 = 0
Theta4 = Desplaza[6, 0]

U5 = 0
V5 = 0
Theta5 = 0

# %% Creación de matrices de desplazamientos y cálculo de los desplazamientos nodales en coordenadas locales


def Vect_des_GLO_and_LOC(ui, vi, thetai, uj, vj, thetaj, Mtrans_Elem):
    """
    Creación de los vectores de desplazamiento para cada uno de los elementos 
    en los cuales se discretiza desplazamientos nodales en el sistema coordenado global como en el sistema coordenado local de cada elemento  
    """
    DesNod_GLO_E = sy.zeros(6, 1)

    DesNod_GLO_E[0, 0] = ui
    DesNod_GLO_E[1, 0] = vi
    DesNod_GLO_E[2, 0] = thetai
    DesNod_GLO_E[3, 0] = uj
    DesNod_GLO_E[4, 0] = vj
    DesNod_GLO_E[5, 0] = thetaj

    DesNod_LOC_E = Mtrans_Elem@DesNod_GLO_E

    return DesNod_GLO_E, DesNod_LOC_E


DesNod_GLO_A, DesNod_LOC_A = Vect_des_GLO_and_LOC(
    U1, V1, Theta1, U2, V2, Theta2, MAtrans)

DesNod_GLO_B, DesNod_LOC_B = Vect_des_GLO_and_LOC(
    U2, V2, Theta2, U3, V3, Theta3, MBtrans)

DesNod_GLO_C, DesNod_LOC_C = Vect_des_GLO_and_LOC(
    U4, V4, Theta4, U3, V3, Theta3, MCtrans)

DesNod_GLO_D, DesNod_LOC_D = Vect_des_GLO_and_LOC(
    U3, V3, Theta3, U5, V5, Theta5, MDtrans)


def Des_locales(DesNod_LOC_Elem):
    """ Creación de las variables de desplazamiento para cada par cada elemento 
    con i,j como Inicio y final del elemento respectivamente. """
    ui_E = DesNod_LOC_Elem[0, 0]
    vi_E = DesNod_LOC_Elem[1, 0]
    thetai_E = DesNod_LOC_Elem[2, 0]
    uj_E = DesNod_LOC_Elem[3, 0]
    vj_E = DesNod_LOC_Elem[4, 0]
    thetaj_E = DesNod_LOC_Elem[5, 0]

    return ui_E, vi_E, thetai_E, uj_E, vj_E, thetaj_E


u1_A, v1_A, theta1_A, u2_A, v2_A, theta2_A = Des_locales(DesNod_LOC_A)

u2_B, v2_B, theta2_B, u3_B, v3_B, theta3_B = Des_locales(DesNod_LOC_B)

u4_C, v4_C, theta4_C, u3_C, v3_C, theta3_C = Des_locales(DesNod_LOC_C)

u3_D, v3_D, theta3_D, u5_D, v5_D, theta5_D = Des_locales(DesNod_LOC_D)

# %% Cálculo de las reacciones

''' 
Reacciones en los nodos 1, 4, 5
'''
"""
REVISAR!!!!
"""
FX1 = float(
    np.array(KA[0, :]@DesNod_GLO_A).astype(np.float64)) + float(VA_Emp_Glo[0, 0])
FY1 = float(
    np.array(KA[1, :]@DesNod_GLO_A).astype(np.float64)) + float(VA_Emp_Glo[1, 0])
M1 = float(np.array(KA[2, :]@DesNod_GLO_A).astype(np.float64)
           ) + float(VA_Emp_Glo[2, 0])

FX4 = float(
    np.array(KC[0, :]@DesNod_GLO_C).astype(np.float64)) + float(VC_Emp_Glo[0, 0])
FY4 = float(
    np.array(KC[1, :]@DesNod_GLO_C).astype(np.float64)) + float(VC_Emp_Glo[1, 0])

FX5 = float(
    np.array(KD[3, :]@DesNod_GLO_D).astype(np.float64)) + float(VD_Emp_Glo[3, 0])
FY5 = float(
    np.array(KD[4, :]@DesNod_GLO_D).astype(np.float64)) + float(VD_Emp_Glo[4, 0])
M5 = float(np.array(KD[5, :]@DesNod_GLO_D).astype(np.float64)
           ) + float(VD_Emp_Glo[5, 0])

dict_reacciones = {"FX1": FX1, "FY1": FY1, "M1": M1,
                   "FX4": FX4, "FY4": FY4, "FX5": FX5, "FY5": FY5, "M5": M5}

Data_reacciones = pd.DataFrame([dict_reacciones]).T

# %% Campos de desplazamientos para cada uno de los elementos

xi = sy.symbols('xi')
# para los elementos tipo portico B, C y D

UH_B, VH_B = cdeshomogeneo(u2_B, u3_B, v2_B, v3_B, theta2_B, theta3_B, x, LB)[
    0], cdeshomogeneo(u2_B, u3_B, v2_B, v3_B, theta2_B, theta3_B, x, LB)[1]
UF_B, VF_B = 0, 0


# El elemento C tiene dos campos de desplazamiento
#   v1 entre 0 < x < 3
#   v2 entre 3 < x < 6
UH_C, VH_C = cdeshomogeneo(u4_C, u3_C, v4_C, v3_C, theta4_C, theta3_C, x, LC)[
    0], cdeshomogeneo(u4_C, u3_C, v4_C, v3_C, theta4_C, theta3_C, x, LC)[1]


#   v1 entre 0 < x < 3

Uf_C1 = sy.integrate(0*Gxx(x, LC, E, Area_portico)[1], (xi, 0, x)) + sy.integrate(
    p_C.subs({x: xi})*Gxx(x, LC, E, Area_portico)[0], (xi, sy.sqrt(10), LC))

Vf_C1 = sy.integrate(0*Gyy(x, xi, LC, E, Ix_P)[1], (xi, 0, x)) + sy.integrate(
    q_C.subs({x: xi})*Gyy(x, xi, LC, E, Ix_P)[0], (xi, sy.sqrt(10), LC))

#   v2 entre 3 < x < 6

Uf_C2 = sy.integrate(p_C.subs({xc: xi})*Gxx(x, LC, E, Area_portico)[1], (xi, sy.sqrt(10), x)) + sy.integrate(
    p_C.subs({xc: xi})*Gxx(x, LC, E, Area_portico)[0], (xi, x, LC))

Vf_C2 = sy.integrate(q_C.subs({xc: xi})*Gyy(x, xi, LC, E, Ix_P)[1], (xi, sy.sqrt(10), x)) + sy.integrate(
    q_C.subs({xc: xi})*Gyy(x, xi, LC, E, Ix_P)[0], (xi, x, LC))

# print('verificación \n', 'derivada 4 = ', sy.diff(sy.simplify(sy.N(Uf_C2)),(x,4))))
# print('Uf_C1 = ', sy.simplify(sy.N(Uf_C1,3)), '\n','Vf_C1 = ', sy.simplify(sy.N(Vf_C1,3)))

UH_D, VH_D = cdeshomogeneo(u3_D, u5_D, v3_D, v5_D, theta3_D, theta5_D, x, LD)[
    0], cdeshomogeneo(u3_D, u5_D, v3_D, v5_D, theta3_D, theta5_D, x, LD)[1]

UF_D, VF_D = cdesempotrado(xi, 0, x, LD, E, Ix_P, Area_portico, p_D.subs({xd: x}), q_D.subs({xd: x}))[
    0], cdesempotrado(xi, 0, x, LD, E, Ix_P, Area_portico, p_D.subs({xd: x}), q_D.subs({xd: x}))[1]

Uf_D = sy.integrate(p_D.subs({xd: xi})*Gxx(x, LD, E, Area_portico)[1], (xi, 0, x)) + sy.integrate(
    (p_D).subs({xd: xi})*Gxx(x, LB, E, Area_portico)[0], (xi, x, LD))

Vf_D = sy.integrate(q_D.subs({xd: xi})*Gyy(x, xi, LD, E, Ix_P)[1], (xi, 0, x)) + sy.integrate(
    (q_D).subs({xd: xi})*Gyy(x, xi, LD, E, Ix_P)[0], (xi, x, LD))
