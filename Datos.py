
import sympy as sy
from Trabajo_estructural_2021 import *
from Graficas_y_Tablas import Elemento_A, Elemento_B, Elemento_C, Elemento_D


with open('Datos_latex.txt', 'w', encoding='utf-8') as f:
    f.write('A continuación se exportan los datos obtenidos mediante el cálculo del método de rigidez')
    f.write('\n\n')

    f.write('\nTransformación de cargas externas de globales a locales\n')
    f.write('Cargas Distribuidas en globales \n')
    f.write('\nElemento C\nQ_C\n')
    f.write(sy.latex(sy.nsimplify(Qc, tolerance=1e-4)))
    f.write('\nElemento D\nQ_C\n')
    f.write(sy.latex(sy.nsimplify(Qd, tolerance=1e-4)))
    f.write('\n\nCargas Distribuidas en locales \n')
    f.write('\nElemento C\np_C o carga axial\n')
    f.write(sy.latex(sy.nsimplify(p_C, tolerance=1e-4)))
    f.write('\nq_C o carga cortante\n')
    f.write(sy.latex(sy.nsimplify(q_C, tolerance=1e-4)))
    f.write('\nElemento D\np_D o carga axial\n')
    f.write(sy.latex(sy.nsimplify(p_D, tolerance=1e-4)))
    f.write('\nq_D o carga cortante\n')
    f.write(sy.latex(sy.nsimplify(q_D, tolerance=1e-4)))

    # Elemento A
    f.write('\n\n\nDatos elemento A\n')
    f.write('Matriz de Rigidez en locales\n')
    f.write((sy.latex(sy.N(MA_loc, 5))))
    f.write('\n\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_A, 5)))
    f.write('\n\nMatriz de transformación\n')
    f.write(sy.latex(sy.N(MAtrans, 5)))
    f.write('\n\nMatriz de Rigidez en coordenadas Globales\n')
    f.write(sy.latex(sy.N(KA, 5)))
    f.write('\n\nVector de empotramiento en globales\n')
    f.write(sy.latex(sy.N(VA_Emp_Glo, 5)))
    f.write('\n\n')

    # Elemento B
    f.write('Datos elemento B\n')
    f.write('Matriz de Rigidez en locales\n')
    f.write((sy.latex(sy.N(MB_loc, 5))))
    f.write('\n\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_B, 5)))
    f.write('\n\nMatriz de transformación\n')
    f.write(sy.latex(sy.N(MBtrans, 5)))
    f.write('\n\nMatriz de Rigidez en coordenadas Globales\n')
    f.write(sy.latex(sy.N(KB, 5)))
    f.write('\n\nVector de empotramiento en globales\n')
    f.write(sy.latex(sy.N(VB_Emp_Glo, 5)))
    f.write('\n\n')

    # Elemento C
    f.write('Datos elemento C\n')
    f.write('Matriz de Rigidez en locales\n')
    f.write((sy.latex(sy.N(MC_loc, 5))))
    f.write('\n\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_C, 5)))
    f.write('\n\nMatriz de transformación\n')
    f.write(sy.latex(sy.N(MCtrans, 5)))
    f.write('\n\nMatriz de Rigidez en coordenadas Globales\n')
    f.write(sy.latex(sy.N(KC, 5)))
    f.write('\n\nVector de empotramiento en globales\n')
    f.write(sy.latex(sy.N(VC_Emp_Glo, 5)))
    f.write('\n\n')

    # Elemento D
    f.write('Datos elemento D\n')
    f.write('Matriz de Rigidez en locales\n')
    f.write((sy.latex(sy.N(MD_loc, 5))))
    f.write('\n\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_D, 5)))
    f.write('\n\nMatriz de transformación\n')
    f.write(sy.latex(sy.N(MDtrans, 5)))
    f.write('\n\nMatriz de Rigidez en coordenadas Globales\n')
    f.write(sy.latex(sy.N(KD, 5)))
    f.write('\n\nVector de empotramiento en globales\n')
    f.write(sy.latex(sy.N(VD_Emp_Glo, 5)))
    f.write('\n\n')

    f.write('Calculo de los desplazamientos deconocidos\n')
    f.write('Sistema de ecuaciones\n')
    f.write(sy.latex(sy.N(Desnod, 5)))
    f.write('\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vect_emp, 5)))
    f.write('\nSolución del sistema de ecuaciones\n\n')
    f.write(sy.latex(des_notacion))
    f.write('\n\n Calculo de las reacciones\n')
    f.write((Data_reacciones.to_latex()))

    # Campos de desplazamiento
    f.write('\n\nCampos de desplazamiento\n')

    # Elemento A
    f.write('Campos de desplazamiento Homogeneo A\n\n')
    f.write('{u^h}_A\n')
    f.write(sy.latex(sy.nsimplify(uh_A, tolerance=1e-5)))
    f.write('\n\n{v^h}_A\n')
    f.write(sy.latex(sy.nsimplify(vh_A, tolerance=1e-5)))
    f.write('\n\nCampos de desplazamientos totales (suma homogeneo más empotrado)\n')
    f.write('u_A\n')
    f.write(sy.latex(sy.nsimplify(uA, tolerance=1e-5)))
    f.write('\n\nv_A\n')
    f.write(sy.latex(sy.nsimplify(vA, tolerance=1e-5)))

    # Elemento B
    f.write('Campos de desplazamiento Homogeneo B\n\n')
    f.write('{u^h}_B\n')
    f.write(sy.latex(sy.nsimplify(uh_B, tolerance=1e-5)))
    f.write('\n\n{v^h}_B\n')
    f.write(sy.latex(sy.nsimplify(vh_B, tolerance=1e-5)))
    f.write('\n\nCampos de desplazamientos totales (suma homogeneo más empotrado)\n')
    f.write('u_B\n')
    f.write(sy.latex(sy.nsimplify(uB, tolerance=1e-5)))
    f.write('\n\nv_B\n')
    f.write(sy.latex(sy.nsimplify(vB, tolerance=1e-5)))

    # Elemento C
    f.write('Campos de desplazamiento Homogeneo C\n\n')
    f.write('{u^h}_C\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uh_C, tolerance=1e-5))))
    f.write('\n\n{v^h}_C\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vh_C, tolerance=1e-5))))

    f.write('\n\n El elemento c posee dos campos de desplazamiento emportrados, ya que la carga externa no se encuentra distribuida en todo el elemento\n')

    # primer tramo
    f.write('\nPrimer tramo\n')
    f.write('{{u^f}_C}^I\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uf_C1, tolerance=1e-5))))
    f.write('\n\n{{v^f}_C}^I\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vf_C1, tolerance=1e-5))))
    f.write('\n\nCampos de desplazamientos totales (suma homogeneo más empotrado)\n')
    f.write('{u_C}^I\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uC_1, tolerance=1e-5))))
    f.write('\n\n{v_C}^I\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vC_1, tolerance=1e-5))))

    # Segundo tramo
    f.write('\n\nSegundo tramo\n')
    f.write('{{u^f}_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uf_C2, tolerance=1e-5))))
    f.write('\n\n{{v^f}_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vf_C2, tolerance=1e-5))))
    f.write('\n\nCampos de desplazamientos totales (suma homogeneo más empotrado)\n')
    f.write('{u_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uC_2, tolerance=1e-5))))
    f.write('\n\n{v_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vC_2, tolerance=1e-5))))

    # Elemento D
    f.write('\n\nCampos de desplazamiento Homogeneo D\n\n')
    f.write('{u^h}_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uh_D, tolerance=1e-5))))
    f.write('\n\n{v^h}_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vh_D, tolerance=1e-5))))
    f.write('\n\nCampos de desplazamiento Empotrado D\n\n')
    f.write('{u^f}_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uf_D, tolerance=1e-5))))
    f.write('\n\n{v^f}_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vf_D, tolerance=1e-5))))
    f.write('\n\nCampos de desplazamientos totales (suma homogeneo más empotrado)\n')
    f.write('u_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(uD, tolerance=1e-5))))
    f.write('\n\nv_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(vD, tolerance=1e-5))))

    # Fuerzas internas
    f.write('\n\n\nFuerzas Internas\n')

    # Elemento A
    f.write('\n\n')
    f.write('Elemento A\n')
    f.write('Fuerza axial\n P_A\n')
    f.write(sy.latex(sy.nsimplify(P_A, tolerance=1e-5)))
    f.write('\nMomento flector\n M_A\n')
    f.write(sy.latex(sy.N(sy.nsimplify(M_A, tolerance=1e-5), 5)))
    f.write('\nFuerza cortante\n V_A\n')
    f.write(sy.latex(sy.N(sy.nsimplify(V_A, tolerance=1e-5), 5)))

    # Elemento B
    f.write('\n\n')
    f.write('Elemento B\n')
    f.write('Fuerza axial\n P_B\n')
    f.write(sy.latex(sy.nsimplify(P_B, tolerance=1e-5)))
    f.write('\nMomento flector\n M_B\n')
    f.write(sy.latex(sy.nsimplify(M_B, tolerance=1e-5)))
    f.write('\nFuerza cortante\n V_B\n')
    f.write(sy.latex(sy.nsimplify(V_B, tolerance=1e-5)))

    # Elemento C
    f.write('\n\n')
    f.write('Elemento C\n')
    f.write('\nPrimer tramo\n')
    f.write('Fuerza axial\n {P_C}^I\n')
    f.write(sy.latex(sy.nsimplify(P_C1, tolerance=1e-5)))
    f.write('\nMomento flector\n {M_C}^I\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(M_C1, tolerance=1e-5))))
    f.write('\nFuerza cortante\n {V_C}^I\n')
    f.write(sy.latex(sy.nsimplify(V_C1, tolerance=1e-5)))
    f.write('\n\nSegundo tramo\n')
    f.write('\nFuerza axial\n {P_C}^II\n')
    f.write(sy.latex(sy.nsimplify(P_C2, tolerance=1e-5)))
    f.write('\nMomento flector\n {M_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(M_C2, tolerance=1e-5))))
    f.write('\nFuerza cortante\n {V_C}^II\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(V_C2, tolerance=1e-5))))

    # Elemento D
    f.write('\n\n')
    f.write('Elemento D\n')
    f.write('Fuerza axial\n P_D\n')
    f.write(sy.latex(sy.nsimplify(P_D, tolerance=1e-5)))
    f.write('\nMomento flector\n M_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(M_D, tolerance=1e-5))))
    f.write('\nFuerza cortante\n V_D\n')
    f.write(sy.latex(sy.expand(sy.nsimplify(V_D, tolerance=1e-5))))

    # Tablas
    f.write('\n\n')
    f.write('Tablas de resumen\n')

    # Elemento A
    f.write('\nElemento A\n')
    f.write((Elemento_A.to_latex(index=False)))

    # Elemento B
    f.write('\nElemento B\n')
    f.write((Elemento_B.to_latex(index=False)))

    # Elemento C
    f.write('\nElemento C\n')
    f.write((Elemento_C.to_latex(index=False)))

    # Elemento D
    f.write('\nElemento D\n')
    f.write((Elemento_D.to_latex(index=False)))

    # Fuerza que el suelo ejerce sobre la viga
    f.write('\nFuerza que el suelo ejerce sobre la viga\n')
    f.write(sy.latex(sy.N(sy.nsimplify(F_soil, tolerance=1e-5), 5)))

    # Fuerzas internas en los extremos de cada elemento

    f.write('\n\n\nElemento A\n')
    f.write('\nMatriz en coordenadas locales\n')
    f.write(sy.latex(sy.N(MA_loc, 5)))
    f.write('\nVector de desplazamientos\n')
    f.write(sy.latex(sy.N(DesNod_LOC_A, 5)))
    f.write('\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_A, 5)))
    f.write('\nResultados\n')
    f.write(sy.latex(sy.N(f_inter_A, 6)))

    f.write('\n\nElemento B\n')
    f.write('\nMatriz en coordenadas locales\n')
    f.write(sy.latex(sy.N(MB_loc, 5)))
    f.write('\nVector de desplazamientos\n')
    f.write(sy.latex(sy.N(DesNod_LOC_B, 5)))
    f.write('\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_B, 5)))
    f.write('\nResultados\n')
    f.write(sy.latex(sy.N(f_inter_B, 6)))

    f.write('\n\nElemento C\n')
    f.write('\nMatriz en coordenadas locales\n')
    f.write(sy.latex(sy.N(MC_loc, 5)))
    f.write('\nVector de desplazamientos\n')
    f.write(sy.latex(sy.N(DesNod_LOC_C, 5)))
    f.write('\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_C, 5)))
    f.write('\nResultados\n')
    f.write(sy.latex(sy.N(f_inter_C, 6)))

    f.write('\n\nElemento D\n')
    f.write('\nMatriz en coordenadas locales\n')
    f.write(sy.latex(sy.N(MD_loc, 5)))
    f.write('\nVector de desplazamientos\n')
    f.write(sy.latex(sy.N(DesNod_LOC_D, 5)))
    f.write('\nVector de empotramiento\n')
    f.write(sy.latex(sy.N(Vec_Emp_D, 5)))
    f.write('\nResultados\n')
    f.write(sy.latex(sy.N(f_inter_D, 6)))
