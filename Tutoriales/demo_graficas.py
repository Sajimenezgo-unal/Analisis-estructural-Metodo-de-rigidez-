# %% Descripción de las graficas
""" Explicación graficas """

# importando las liberias necesarias
import sympy as sy
import matplotlib.pyplot as plt
from matplotlib import transforms

# Creando los valores par x correspondientes a cada elemento
# Donde LA, LB, LD, LC corresponden a las longitudes de los elementos

xA = np.linspace(0, LA, 100)
xB = np.linspace(0, LB, 100)
xD = np.linspace(0, float(LD), 100)

# Elemento c la carga comienza en sqrt(5)

# xC1 = np.linspace(0,'valor inicio carga',100)
# xC2 = np.linspace('valor inicio carga','longitud de C', 100)

xC1 = np.linspace(0, float(sy.sqrt(5)), 100)
xC2 = np.linspace(float(sy.sqrt(5)), float(LC), 100)

Escala = 100  # Escala ampliada 100 veces
EscalaP = 0.01  # Escala 1/100
EscalaV = 0.001  # Escala 1/1000
EscalaM = 0.01  # Escala 1/1000


# Se crean las variables donde serán almacenados los resultados de los
# campos evaluados posteriormente

# Elemento A
uA_gra, vA_gra = [], []
PA_gra, MA_gra, VA_gra = [], [], []
# Elemento B
uB_gra, vB_gra = [], []
PB_gra, MB_gra, VB_gra = [], [], []
# Elemento C
uC1_gra, vC1_gra = [], []
uC2_gra, vC2_gra = [], []
PC1_gra, PC2_gra = [], []
MC1_gra, MC2_gra = [], []
VC1_gra, VC2_gra = [], []
# Elemento D
uD_gra, vD_gra = [], []
PD_gra, MD_gra, VD_gra = [], [], []

# Evaluación de los campos del elemento A
for i in xA:
    uA_gra.append(u_A.subs({x: i})*Escala)
    vA_gra.append(v_A.subs({x: i})*Escala)
    PA_gra.append(P_A.subs({x: i})*EscalaP)
    MA_gra.append(M_A.subs({x: i})*EscalaM)
    VA_gra.append(V_A.subs({x: i})*EscalaV)

# Evaluación de los campos del elemento B
for i in xB:
    uB_gra.append(u_B.subs({x: i})*Escala)
    vB_gra.append(v_B.subs({x: i})*Escala*Escala)
    PB_gra.append(P_B.subs({x: i})*EscalaP)
    MB_gra.append(M_B.subs({x: i})*EscalaM)
    VB_gra.append(V_B.subs({x: i})*EscalaV)

# Evaluación de los campos del elemento D
for i in xD:
    uD_gra.append(u_D.subs({x: i})*Escala)
    vD_gra.append(v_D.subs({x: i})*Escala)
    PD_gra.append(P_D.subs({x: i})*EscalaP)
    MD_gra.append(M_D.subs({x: i})*EscalaM)
    VD_gra.append(V_D.subs({x: i})*EscalaV)

# Evaluación de los campos del elemento C
#   Antes de la carga
for i in xC1:
    uC1_gra.append(u_C1.subs({x: i})*Escala)
    vC1_gra.append(v_C1.subs({x: i})*Escala)
    PC1_gra.append(P_C1.subs({x: i})*EscalaP)
    MC1_gra.append(M_C1.subs({x: i})*EscalaM)
    VC1_gra.append(V_C1.subs({x: i})*EscalaV)
    ##############################
#   Despues de la carga
for i in xC2:
    uC2_gra.append(u_C2.subs({x: i})*Escala)
    vC2_gra.append(v_C2.subs({x: i})*Escala)
    PC2_gra.append(P_C2.subs({x: i})*EscalaP)
    MC2_gra.append(M_C2.subs({x: i})*EscalaM)
    VC2_gra.append(V_C2.subs({x: i})*EscalaV)

x_elem_C = np.linspace(0, 6, 100)

# Elem_C = []
# for i in x_elem_C:
#     y = (2/6)*x + 7
#     Elem_C.append(y.subs({x: i}))

# x_elem_D = np.linspace(0, 5, 100)

# Elem_D = []
# for i in x_elem_D:
#     y = -(2/5)*x + 9
#     Elem_D.append(y.subs({x: i}))

# Transformación de los campos de desplazamiento y fuerzas internas
#   rot_C = transforms.Affine2D().rotate('angulo en radianes').translate(x,y)
#   donde (x,y) punto de incio de la grafica rotada

base = plt.gca().transData
rot_C = transforms.Affine2D().rotate(theta_C).translate(0, 7)
rot_A = transforms.Affine2D().rotate(radians(90)).translate(6, 0)
rot_B = transforms.Affine2D().rotate(radians(90)).translate(6, 5)
rot_D = transforms.Affine2D().rotate(theta_D).translate(6, 9)

# Grafica

# Formato de LaTex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Grafica de la estructura
# plt.axvline(x, yinicio, yfinal/10, color='negro(k)',linewidth= 1(grosor de la linea))
plt.axvline(6, 0, 0.9, color='k', linewidth=1)  # linea Elemento A y B
plt.axhline(y=5, color='k', linestyle=':', linewidth=0.8)  # linea suelo en y=5


# grafica campos de desplazamiento
plt.plot(xA, uA_gra, color='k', linestyle='--', transform=rot_A + base)
plt.plot(xB, uB_gra, color='k', linestyle='--', transform=rot_B + base)
plt.plot(xC1, uC1_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xC2, uC2_gra, color='k', linestyle='--', transform=rot_C + base)
plt.plot(xD, uD_gra, color='k', linestyle='--', transform=rot_D + base)

plt.axis([-1, 12, 0, 10])
plt.xlabel('$x$[m]', fontsize=16, weight="bold")
plt.ylabel('$y$[m]', fontsize=16, weight="bold")
plt.grid(linestyle='--', linewidth=0.3, color='k')
plt.show()
# plt.savefig('Campo de desplazamiento horizontal.pdf', bbox_inches='tight')
# %%
