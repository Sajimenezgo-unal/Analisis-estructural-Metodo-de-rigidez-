
"""
Dentro del siguiente codigo se realizan las graficas y tablas requeridas para el trabajo 
estas se separan del resto del código con el fin de llevar un desarrollo más ordenado 
en archivos menos densos para el lector 
"""

from Modulos import *
from Trabajo_estructural_2021 import *

####
# Librerias
import sympy as sy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import transforms
from sympy import init_printing
init_printing()

# %% Graficas

# LA = 5  # longitud de A en metros
# LB = 4  # longitud de B en metros
# LC = sy.sqrt(6**2 + 2**2)  # longitud de C en metros
# LD = sy.sqrt(5**2 + 2**2)  # longitud de D en metros


xA_eval = np.linspace(0, LA, 100)
xB_eval = np.linspace(0, LB, 100)
xC_eval1 = np.linspace(0, float(sy.sqrt(10)), 100)
xC_eval2 = np.linspace(float(sy.sqrt(10)), float(LC), 100)
xD_eval = np.linspace(0, float(LD), 100)

Escala = 100
EscalaP = 0.01
EscalaV = 0.001
EscalaM = 0.01

# Datos Graficas

# Elemento A
u_A_gra, v_A_gra = [], []
P_A_gra = []
M_A_gra = []
V_A_gra = []
# Elemento B
u_B_gra, v_B_gra = [], []
P_B_gra = []
M_B_gra = []
V_B_gra = []
# Elemento C
u_C1_gra, v_C1_gra = [], []
u_C2_gra, v_C2_gra = [], []
P_C1_gra, P_C2_gra = [], []
M_C1_gra, M_C2_gra = [], []
V_C1_gra, V_C2_gra = [], []
# Elemento D
u_D_gra, v_D_gra = [], []
P_D_gra = []
M_D_gra = []
V_D_gra = []

# Elemento A

for i in xA_eval:
    u_A_gra.append(uA.subs({x: i})*Escala)
    v_A_gra.append(vA.subs({x: i})*Escala)
    P_A_gra.append(P_A.subs({x: i})*EscalaP)
    M_A_gra.append(M_A.subs({x: i})*EscalaM)
    V_A_gra.append(V_A.subs({x: i})*EscalaV)
# Elemento B

for i in xB_eval:
    u_B_gra.append(uB.subs({x: i})*Escala)
    v_B_gra.append(vB.subs({x: i})*Escala)
    P_B_gra.append(P_B.subs({x: i})*EscalaP)
    M_B_gra.append(M_B.subs({x: i})*EscalaM)
    V_B_gra.append(V_B.subs({x: i})*EscalaV)

# Elemento C

for i in xC_eval1:
    u_C1_gra.append(uC_1.subs({x: i})*Escala)
    v_C1_gra.append(vC_1.subs({x: i})*Escala)
    P_C1_gra.append(P_C1.subs({x: i})*EscalaP)
    M_C1_gra.append(M_C1.subs({x: i})*EscalaM)
    V_C1_gra.append(V_C1.subs({x: i})*EscalaV)
    ###########################
for i in xC_eval2:
    u_C2_gra.append(uC_2.subs({x: i})*Escala)
    v_C2_gra.append(vC_2.subs({x: i})*Escala)
    P_C2_gra.append(P_C2.subs({x: i})*EscalaP)
    M_C2_gra.append(M_C2.subs({x: i})*EscalaM)
    V_C2_gra.append(V_C2.subs({x: i})*EscalaV)

# Elemento D

for i in xD_eval:
    u_D_gra.append(uD.subs({x: i})*Escala)
    v_D_gra.append(vD.subs({x: i})*Escala)
    P_D_gra.append(P_D.subs({x: i})*EscalaP)
    M_D_gra.append(M_D.subs({x: i})*EscalaM)
    V_D_gra.append(V_D.subs({x: i})*EscalaV)

x_elem_C = np.linspace(0, 6, 100)

Elem_C = []
for i in x_elem_C:
    y = (2/6)*x + 7
    Elem_C.append(y.subs({x: i}))

x_elem_D = np.linspace(0, 5, 100)

Elem_D = []
for i in x_elem_D:
    y = -(2/5)*x + 9
    Elem_D.append(y.subs({x: i}))

"""
Transformación de la grafica utilizando 
base = plt.gca().transData
rot = transforms.Affine2D().rotate('angle in radians').translate(x,y)
"""

base = plt.gca().transData
rot_C = transforms.Affine2D().rotate(theta_C).translate(0, 7)
rot_A = transforms.Affine2D().rotate(radians(90)).translate(6, 0)
rot_B = transforms.Affine2D().rotate(radians(90)).translate(6, 5)
rot_D = transforms.Affine2D().rotate(theta_D).translate(6, 9)

# %%              campo de desplazamiento vertical

# Estructura

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.axvline(6, 0, 0.9, color='k', linewidth=0.8)
plt.axhline(y=5, color='k', label='Suelo', linestyle=':', linewidth=0.8)
plt.plot(x_elem_C, Elem_C, color='k', linewidth=0.8)
plt.plot(np.linspace(6, 11, 100), Elem_D, color='k', linewidth=0.8)

# grid
plt.axvline(1, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=1, color='k', linestyle='-', linewidth=0.3)
plt.axvline(3, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=3, color='k', linestyle='--', linewidth=0.3)
plt.axvline(5, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axvline(7, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=7, color='k', linestyle='--', linewidth=0.3)
plt.axvline(9, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=9, color='k', linestyle='--', linewidth=0.3)
plt.axvline(11, 0, 11, color='k', linestyle='--', linewidth=0.3)

# Campos
plt.plot(xC_eval1, v_C1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC_eval2, v_C2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xA_eval, v_A_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB_eval, v_B_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xD_eval, v_D_gra, color='k', linestyle='--', transform=rot_D + base)

# # values
# plt.text(np.array(7), np.array(0.5), f"{sy.N(v_A_gra[0],3)}")
# plt.text(np.array(7), np.array(4.5), f"{sy.N(v_A_gra[99],3)}")

# plt.text(np.array(7), np.array(5.5), f"{sy.N(v_B_gra[0],3)}")
# plt.text(np.array(6), np.array(9.5), f"{sy.N(v_B_gra[99],3)}")

# plt.text(np.array(-0.3), np.array(7.5), f"{sy.N(v_C1_gra[0],3)}")
# plt.text(np.array(4.2), np.array(9.3), f"{sy.N(v_C2_gra[99],3)}")

# plt.text(np.array(7.3), np.array(8.7), f"{sy.N(v_D_gra[0],3)}")
# plt.text(np.array(10), np.array(6), f"{sy.N(v_D_gra[99],3)}")


plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16)
plt.ylabel('$y$[m]', fontsize=16)
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.savefig('LaTex/Graficas/Grafica_Des_vertical.pdf', bbox_inches='tight')
# plt.show()
plt.close()
# %%             campo de desplazamiento horizontal

# Estructura

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.axvline(6, 0, 0.5, label='Elemento A', color='k', linewidth=0.8)
plt.axvline(6, 0.5, 0.9, label='Elemento B', color='k', linewidth=0.8)
plt.axhline(y=5, color='k', linestyle=':', linewidth=0.8)
plt.plot(x_elem_C, Elem_C, label='Elemento C', color='k', linewidth=0.8)
plt.plot(np.linspace(6, 11, 100), Elem_D,
         label='Elemento D', color='k', linewidth=0.8)

# grid
plt.axvline(1, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=1, color='k', linestyle='-', linewidth=0.3)
plt.axvline(3, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=3, color='k', linestyle='--', linewidth=0.3)
plt.axvline(5, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axvline(7, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=7, color='k', linestyle='--', linewidth=0.3)
plt.axvline(9, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=9, color='k', linestyle='--', linewidth=0.3)
plt.axvline(11, 0, 11, color='k', linestyle='--', linewidth=0.3)

# Campos
plt.plot(xC_eval1, u_C1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC_eval2, u_C2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xA_eval, u_A_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB_eval, u_B_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xD_eval, u_D_gra, color='k', linestyle='--', transform=rot_D + base)

# # values
# plt.text(np.array(xA_eval[0] + 0.3), np.array(u_A_gra[0] + 0.3), f"{sy.N(u_A_gra[0],3)}",
#                                                                  transform=rot_A + base)
# plt.text(np.array(xA_eval[99] + 0.3), np.array(u_A_gra[99] + 0.3), f"{sy.N(u_A_gra[99],3)}",
#                                                                  transform=rot_A + base)

# plt.text(np.array(xB_eval[0] + 0.3), np.array(u_B_gra[0] + 0.3), f"{sy.N(u_B_gra[0],3)}",
#                                                                  transform=rot_B + base)
# plt.text(np.array(xB_eval[99] + 0.3), np.array(u_B_gra[99] + 0.3), f"{sy.N(u_B_gra[99],3)}",
#                                                                  transform=rot_B + base)

# plt.text(np.array(xC_eval1[0] + 0.3), np.array(u_C1_gra[0] + 0.3), f"{sy.N(u_C1_gra[0],3)}",
#                                                                  transform=rot_C + base)
# plt.text(np.array(xC_eval2[99] + 0.3), np.array(u_C2_gra[99] + 0.3), f"{sy.N(u_C2_gra[99],3)}",
#                                                                  transform=rot_C + base)

# plt.text(np.array(xD_eval[0] + 0.3), np.array(u_D_gra[0] + 0.3), f"{sy.N(u_D_gra[0],3)}",
#                                                                  transform=rot_D + base)
# plt.text(np.array(xD_eval[99] + 0.3), np.array(u_D_gra[99] + 0.3), f"{sy.N(u_D_gra[99],3)}",
#                                                                  transform=rot_D + base)


plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16, weight="bold")
plt.ylabel('$y$[m]', fontsize=16, weight="bold")
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.savefig('LaTex/Graficas/Grafica_Des_horizontal.pdf', bbox_inches='tight')
# plt.show()
plt.close()

# %%              Campo de fuerza axial

# Estructura

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.axvline(6, 0, 0.5, label='Elemento A', color='k', linewidth=0.8)
plt.axvline(6, 0.5, 0.9, label='Elemento B', color='k', linewidth=0.8)
plt.axhline(y=5, color='k', linestyle=':', linewidth=0.8)
plt.plot(x_elem_C, Elem_C, label='Elemento C', color='k', linewidth=0.8)
plt.plot(np.linspace(6, 11, 100), Elem_D,
         label='Elemento D', color='k', linewidth=0.8)

# grid
plt.axvline(1, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=1, color='k', linestyle='-', linewidth=0.3)
plt.axvline(3, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=3, color='k', linestyle='--', linewidth=0.3)
plt.axvline(5, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axvline(7, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=7, color='k', linestyle='--', linewidth=0.3)
plt.axvline(9, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=9, color='k', linestyle='--', linewidth=0.3)
plt.axvline(11, 0, 11, color='k', linestyle='--', linewidth=0.3)

# Campos
plt.plot(xC_eval1, P_C1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC_eval2, P_C2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xA_eval, P_A_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB_eval, P_B_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xD_eval, P_D_gra, color='k', linestyle='--', transform=rot_D + base)

# # values
plt.text(np.array(xA_eval[0] + 0.5), np.array(P_A_gra[0] - 0.5),
         f"{sy.N(P_A_gra[0],3)}", transform=rot_A + base)
# plt.text(np.array(xA_eval[99] + 0.3), np.array(P_A_gra[99] + 0.3),
#          f"{sy.nsimplify(P_A_gra[99],tolerance=1e-6)}", transform=rot_A + base)

# plt.text(np.array(xB_eval[0] - 0.5), np.array(P_B_gra[0] + 0.3),
#          f"{sy.N(P_B_gra[0],3)}", transform=rot_B + base)
plt.text(np.array(xB_eval[99] + 0.3), np.array(P_B_gra[99] +
                                               0.3), f"{sy.N(P_B_gra[99],3)}", transform=rot_B + base)

plt.text(np.array(xC_eval1[0] - 1), np.array(P_C1_gra[0] + 0.3),
         f"{sy.N(P_C1_gra[0],3)}", transform=rot_C + base)
plt.text(np.array(xC_eval2[99] - 2), np.array(P_C2_gra[99] + 0.4),
         f"${sy.N(P_C2_gra[99],3)}^C$", transform=rot_C + base)

plt.text(np.array(xD_eval[0] + 0.3), np.array(P_D_gra[0] + 0.5),
         f"${sy.N(P_D_gra[0],3)}^D$", transform=rot_D + base)
plt.text(np.array(xD_eval[99] - 0.3), np.array(P_D_gra[99] -
                                               0.5), f"{sy.N(P_D_gra[99],3)}", transform=rot_D + base)


plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16, weight="bold")
plt.ylabel('$y$[m]', fontsize=16, weight="bold")
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.savefig('LaTex/Graficas/Grafica_campo_fuerza_axial.pdf',
            bbox_inches='tight')
# plt.show()
plt.close()

# %%              Campo de momento flector

# Estructura

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.axvline(6, 0, 0.5, label='Elemento A', color='k', linewidth=0.8)
plt.axvline(6, 0.5, 0.9, label='Elemento B', color='k', linewidth=0.8)
plt.axhline(y=5, color='k', linestyle=':', linewidth=0.8)
plt.plot(x_elem_C, Elem_C, label='Elemento C', color='k', linewidth=0.8)
plt.plot(np.linspace(6, 11, 100), Elem_D,
         label='Elemento D', color='k', linewidth=0.8)

# grid
plt.axvline(1, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=1, color='k', linestyle='-', linewidth=0.3)
plt.axvline(3, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=3, color='k', linestyle='--', linewidth=0.3)
plt.axvline(5, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axvline(7, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=7, color='k', linestyle='--', linewidth=0.3)
plt.axvline(9, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=9, color='k', linestyle='--', linewidth=0.3)
plt.axvline(11, 0, 11, color='k', linestyle='--', linewidth=0.3)

# Campos
plt.plot(xC_eval1, M_C1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC_eval2, M_C2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xA_eval, M_A_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB_eval, M_B_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xD_eval, M_D_gra, color='k', linestyle='--', transform=rot_D + base)

# # values
plt.text(np.array(xA_eval[0] + 0.5), np.array(M_A_gra[0] - 0.3),
         f"{sy.N(M_A_gra[0],3)}", transform=rot_A + base)
# plt.text(np.array(xA_eval[99] - 0.5), np.array(M_A_gra[99] - 0.5), f"{sy.N(M_A_gra[99],3)}",transform=rot_A + base)

plt.text(np.array(xB_eval[0] + 0.3), np.array(M_B_gra[0] - 0.5),
         f"{sy.N(M_B_gra[0],3)}", transform=rot_B + base)
plt.text(np.array(xB_eval[99] + 0.5), np.array(M_B_gra[99] - 0.3),
         f"${sy.N(M_B_gra[99],3)}^B$", transform=rot_B + base)

plt.text(np.array(xC_eval1[0] - 1), np.array(M_C1_gra[0] + 0.3),
         f"{sy.N(M_C1_gra[0],3)}", transform=rot_C + base)
plt.text(np.array(xC_eval2[99] - 0.3), np.array(M_C2_gra[99] - 0.5),
         f"${sy.N(M_C2_gra[99],3)}^C$", transform=rot_C + base)

plt.text(np.array(xD_eval[0] - 1.5), np.array(M_D_gra[0]-1),
         f"${sy.N(M_D_gra[0],3)}^D$", transform=rot_D + base)
plt.text(np.array(xD_eval[99] - 0.3), np.array(M_D_gra[99] -
                                               0.5), f"{sy.N(M_D_gra[99],3)}", transform=rot_D + base)


plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16, weight="bold")
plt.ylabel('$y$[m]', fontsize=16, weight="bold")
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.savefig('LaTex/Graficas/Grafica_campo_momento_flector.pdf',
            bbox_inches='tight')
# plt.show()
plt.close()

# %%              Campo de fuerza cortante

# Estructura

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.axvline(6, 0, 0.5, label='Elemento A', color='k', linewidth=0.8)
plt.axvline(6, 0.5, 0.9, label='Elemento B', color='k', linewidth=0.8)
plt.axhline(y=5, color='k', linestyle=':', linewidth=0.8)
plt.plot(x_elem_C, Elem_C, label='Elemento C', color='k', linewidth=0.8)
plt.plot(np.linspace(6, 11, 100), Elem_D,
         label='Elemento D', color='k', linewidth=0.8)
# plt.text(6.3, 2.5, 'A')
# plt.text(6.3, 7, 'B')
# plt.text(3.4, 7.3, 'C')
# plt.text(8.5, 7.3, 'D')

# grid
plt.axvline(1, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=1, color='k', linestyle='-', linewidth=0.3)
plt.axvline(3, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=3, color='k', linestyle='--', linewidth=0.3)
plt.axvline(5, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axvline(7, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=7, color='k', linestyle='--', linewidth=0.3)
plt.axvline(9, 0, 11, color='k', linestyle='--', linewidth=0.3)
plt.axhline(y=9, color='k', linestyle='--', linewidth=0.3)
plt.axvline(11, 0, 11, color='k', linestyle='--', linewidth=0.3)

# Campos
plt.plot(xC_eval1, V_C1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC_eval2, V_C2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xA_eval, V_A_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB_eval, V_B_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xD_eval, V_D_gra, color='k', linestyle='--', transform=rot_D + base)


# values
plt.text(np.array(xA_eval[0]+0.5), np.array(V_A_gra[0]-0.5),
         f"{sy.N(V_A_gra[0],3)}", transform=rot_A + base)
# plt.text(np.array(xA_eval[99]-0.5), np.array(V_A_gra[99]-0.5), f"{sy.N(V_A_gra[99],3)}",transform=rot_A + base)

plt.text(np.array(xB_eval[0]), np.array(V_B_gra[0]-0.5),
         f"{sy.N(V_B_gra[0],3)}", transform=rot_B + base)
plt.text(np.array(xB_eval[99]+0.4), np.array(V_B_gra[99]-0.3),
         f"${sy.N(V_B_gra[99],3)}^B$", transform=rot_B + base)

plt.text(np.array(xC_eval1[0]-0.3), np.array(V_C1_gra[0]-0.7),
         f"{sy.N(V_C1_gra[0],3)}", transform=rot_C + base)
plt.text(np.array(xC_eval2[99]-1.5), np.array(V_C2_gra[99]+0.5),
         f"{sy.N(V_C2_gra[99],3)}", transform=rot_C + base)

plt.text(np.array(xD_eval[0]+0.5), np.array(V_D_gra[0]-1),
         f"${sy.N(V_D_gra[0],3)}^D$", transform=rot_D + base)
plt.text(np.array(xD_eval[99]-0.8), np.array(V_D_gra[99]-0.8),
         f"{sy.N(V_D_gra[99],3)}", transform=rot_D + base)

plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16, weight="bold")
plt.ylabel('$y$[m]', fontsize=16, weight="bold")
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.savefig('LaTex/Graficas/Grafica_campo_fuerza_cortante.pdf',
            bbox_inches='tight')
# plt.show()
plt.close()


# %% Tablas

def Tablas_fuerzas_internas(FI, xi, xj, n, F, Elem):
    """
    Genera las tablas para los valores de x_E/10
    Donde:  
            FI: Fuerza interna
            xi,xj: intervalo de evaluación de FI
            n: número de veces que se divide el intervalo
            F: nombre de la fuerza interna
            Elem: Elemento al que coresponde la tabla 
    """

    import pandas as pd
    xE = np.linspace(float(xi), float(xj), n)
    data = {f"$x_{Elem}$": [], f"${F}$": []}
    for i in xE:
        valor = np.round(float(FI.subs({x: i})), decimals=5)
        data[f"$x_{Elem}$"].append(np.round(float(i), decimals=5))
        data[f"${F}$"].append(valor)
    tabla = pd.DataFrame.from_dict(data)

    return tabla


A, B, C, D = sy.symbols('A,B,C,D')
#    Elemento A
PA_Tabla = Tablas_fuerzas_internas(P_A, 0, LA, 10, "P($x'_A$)", A)
MA_Tabla = Tablas_fuerzas_internas(M_A, 0, LA, 10, "M($x'_A$)", A)
VA_Tabla = Tablas_fuerzas_internas(V_A, 0, LA, 10, "V($x'_A$)", A)
#   Elemento B
PB_Tabla = Tablas_fuerzas_internas(P_B, 0, LB, 10, "P($x'_B$)", B)
MB_Tabla = Tablas_fuerzas_internas(M_B, 0, LB, 10, "M($x'_B$)", B)
VB_Tabla = Tablas_fuerzas_internas(V_B, 0, LB, 10, "V($x'_B$)", B)
#   Elemento D
PD_Tabla = Tablas_fuerzas_internas(P_D, 0, float(LD), 10, "P($x'_C$)", D)
MD_Tabla = Tablas_fuerzas_internas(M_D, 0, float(LD), 10, "M($x'_C$)", D)
VD_Tabla = Tablas_fuerzas_internas(V_D, 0, float(LD), 10, "V($x'_C$)", D)
#   Elemento C
PC1_Tabla = Tablas_fuerzas_internas(P_C1, 0, float(
    LC), 10, "P($x'_C$)", C).drop(index=[5, 6, 7, 8, 9])
MC1_Tabla = Tablas_fuerzas_internas(M_C1, 0, float(
    LC), 10, "M($x'_C$)", C).drop(index=[5, 6, 7, 8, 9])
VC1_Tabla = Tablas_fuerzas_internas(V_C1, 0, float(
    LC), 10, "V($x'_C$)", C).drop(index=[5, 6, 7, 8, 9])
PC2_Tabla = Tablas_fuerzas_internas(P_C2, 0, float(
    LC), 10, "P($x'_C$)", C).drop(index=[0, 1, 2, 3, 4])
MC2_Tabla = Tablas_fuerzas_internas(M_C2, 0, float(
    LC), 10, "M($x'_C$)", C).drop(index=[0, 1, 2, 3, 4])
VC2_Tabla = Tablas_fuerzas_internas(V_C2, 0, float(
    LC), 10, "V($x'_C$)", C).drop(index=[0, 1, 2, 3, 4])

#       Merge fuerzas internas Elemento C

PC_Tabla = pd.concat([PC1_Tabla, PC2_Tabla])
MC_Tabla = pd.concat([MC1_Tabla, MC2_Tabla])
VC_Tabla = pd.concat([VC1_Tabla, VC2_Tabla])

Elemento_A = pd.merge(PA_Tabla, pd.merge(
    VA_Tabla, MA_Tabla, on='$x_A$'), on='$x_A$')
Elemento_B = pd.merge(PB_Tabla, pd.merge(
    VB_Tabla, MB_Tabla, on='$x_B$'), on='$x_B$')
Elemento_C = pd.merge(PC_Tabla, pd.merge(
    VC_Tabla, MC_Tabla, on='$x_C$'), on='$x_C$')
Elemento_D = pd.merge(PD_Tabla, pd.merge(
    VD_Tabla, MD_Tabla, on='$x_D$'), on='$x_D$')

Elemento_A.name = 'Elemento_A'
Elemento_B.name = 'Elemento_B'
Elemento_C.name = 'Elemento_C'
Elemento_D.name = 'Elemento_D'
