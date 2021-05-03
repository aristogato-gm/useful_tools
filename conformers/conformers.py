import MDAnalysis as mda
from itertools import product
from myconfig import * 
import numpy as np

#leer topologia y trajectoria
u = mda.Universe("file.gro","file.trr")
print(u)
print("La trayectoria tiene {} frames".format(len(u.trajectory)))
print("El número de moléculas es {}".format(len(u.atoms.residues)))
num_mol = len(u.atoms.residues)   

print("El número de diedros que se calculará es: {}".format(num_diedros))

diedros_list = []
for j in range(1,num_diedros+1):
    diedros_list.append(vars()['diedro{}'.format(j)])
    print("El diedro {} tiene los siguientes índices {}".format(j,diedros_list[j-1]))

print(diedros_list)
#print(diedros_list[2])

la_lista = []
for i in range(1,num_mol+1):
    for ts in u.trajectory[inicio:fin+1]:
        time = u.trajectory.time
        UNK3 = u.select_atoms('resid {}'.format(i))
#seleccion del diedro de interés
        die_list = []
        for j in range(0,num_diedros):
            die_list.append(UNK3.atoms[diedros_list[j]])
#calcular el ángulo diedro y ajustarlo para que su valor esté entre 0 y 360 grados
            if die_list[j].dihedral.value()<0:
                die_list[j] = 360+die_list[j].dihedral.value()
            else:
                die_list[j]=die_list[j].dihedral.value()
#clasificar el ángulo diedro como t, g o gp
            if die_list[j]>120 and die_list[j]<240:
                die_list[j] = 't'
            elif die_list[j]>0 and die_list[j]<120:
                die_list[j] = 'g'
            elif die_list[j]>240 and die_list[j]<360:
                die_list[j] = 'gp'

        tupla1=tuple(die_list)
        la_lista.append(tupla1)
#print(la_lista)    
my_array = np.asarray(la_lista)
#contar confórmeros creando permutaciones con repetición de num_diedros
print("EL PORCENTAJE DE CONFÓRMEROS ES EL SIGUIENTE:")
perm2 = ['t', 'g', 'gp']
suma = 0
for roll in product(list(perm2), repeat = num_diedros):
    contador = 0
    for i in la_lista:
        if i==roll:
            contador=contador+1
    suma = suma+contador
    porcentaje = (contador*100)/(num_mol*(fin-inicio))
    print("Confórmero {0} hay: {1:.2f}% ".format(roll,porcentaje))
print(suma/(fin-inicio))
