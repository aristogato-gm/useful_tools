import MDAnalysis as mda
import numpy as np
import MDAnalysis.transformations

#decide across what axis we want join our two systems,
#the options are: axis = 0 for x axis, axis = 1 for y and axis = 2 for z)
axis = 2

#gap between the two systems in Angstroms
gap = 10

#load the two structure to create the universe u1 and u2
u1 = mda.Universe("system1.gro")
u2 = mda.Universe("system2.gro")

atoms1 = u1.select_atoms('resname *')
posiciones1 = atoms1.positions
centroid1 = atoms1.centroid()

atoms2 = u2.select_atoms('resname *')
posiciones2 = atoms2.positions
centroid2 = atoms2.centroid()

delta = centroid1 - centroid2
delta[axis] = centroid1[axis]*2 + gap

workflow = [MDAnalysis.transformations.translate(delta)]
u2 = mda.Universe("system2.gro", transformations=workflow)

cell_size = u1.dimensions
cell_size[axis] = u1.dimensions[axis] + u2.dimensions[axis] + gap

u3 = mda.Merge(u1.atoms, u2.atoms)
u3.dimensions = cell_size

u3.atoms.write('out.gro')
