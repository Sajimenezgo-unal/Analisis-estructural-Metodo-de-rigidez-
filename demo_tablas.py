# %% tablas

"""
Comentarios: usar print al exportar las tablas para 
            evitar el error por saltos de linea (\\\\\)
            Revisar los de empalme de las funciones de ser necesarios
            en caso contrario eliminar el .drop(index[])
            Revisar en que punto (indice) de  los valores de x cambia de 
            funcion para soltar los correctos
"""

# se importan las funciones del código Base
from Trabajo_estructural_2021 import *
import pandas as pd

# Se crea la función que desarrolla las tablas


def tablas(FI, x0, x1, n, F, Elem):
    """
    Genera las tablas para los valores de x_E/10
    Donde:  
            FI: Fuerza interna
            x0,x1: intervalo de evaluación de FI
            n: número de veces que se divide el intervalo
            F: nombre de la fuerza interna
            Elem: Elemento al que coresponde la tabla 
    """
    xvalores = np.linspace(float(x0), float(x1), n)  # valores de x a evaluar
    """
    Diccionario = {'nombre1':'objeto1','nombre2':'objeto2'},
    donde los objetos pueden ser valores o listas o diccionario 
    """
    Datos = {f'$x_{Elem}$': [], f'${F}$': []}
    for xn in xvalores:

        valor = np.round(float(FI.subs({x: xn})), decimals=5)
        Datos[f'$x_{Elem}$'].append(np.round(float(xn), decimals=5))
        Datos[f'${F}$'].append(valor)
    tabla = pd.DataFrame.from_dict(Datos)
    return tabla

# Elemento C el cual cuenta con dos funciones por funciones de fuerza interna


# primer tramos c1
# los indices pueden variar
tabla_Pc1 = tablas(P_C1, 0, float(LC), 50, 'P(x_c)[kN]', 'C').drop(
    index=[5, 6, 7, 8, 9])
tabla_Mc1 = tablas(M_C1, 0, float(LC), 50, 'M(x_c)[kN]', 'C').drop(
    index=[5, 6, 7, 8, 9])
tabla_Vc1 = tablas(V_C1, 0, float(LC), 50, 'V(x_c)[kN]', 'C').drop(
    index=[5, 6, 7, 8, 9])
# primer tramos c2
# los indices pueden variar
tabla_Pc2 = tablas(P_C2, 0, float(LC), 50, 'P(x_c)[kN]', 'C').drop(
    index=[0, 1, 2, 3, 4])
tabla_Mc2 = tablas(M_C2, 0, float(LC), 50, 'M(x_c)[kN]', 'C').drop(
    index=[0, 1, 2, 3, 4])
tabla_Vc2 = tablas(V_C2, 0, float(LC), 50, 'V(x_c)[kN]', 'C').drop(
    index=[0, 1, 2, 3, 4])

# Se utiliza el método drop para soltar o eliminar las parte de las tablas con
# funciones que fueron evaluadas pero que no corresponden a las funciones correctas a evaluar

# Se concatenan las dos porciones de la tabla con el fin de tener una tabla con 10 valores de x evaluados
tabla_PC = pd.concat([tabla_Pc1, tabla_Pc2])
tabla_MC = pd.concat([tabla_Mc1, tabla_Mc2])
tabla_VC = pd.concat([tabla_Vc1, tabla_Vc2])

# Merge para juntar las tablas en una sola por elemento
tabla_C = pd.merge(tabla_PC, pd.merge(
    tabla_MC, tabla_VC, on='$x_C$'), on='$x_C$')

# exportar a latex
# nombretabla.to_latex(index=False)
print(tabla_C.to_latex(index=False))

# %%
