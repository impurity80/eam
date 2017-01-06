from ase.lattice.cubic import *
from ase.lattice.surface import surface
from ase.calculators.eam import EAM
from ase.visualize import view, write
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
import elastic
from ase.optimize import BFGS
from ase.build import fcc111

from ase.md.langevin import Langevin
from ase import units
from ase.io.trajectory import Trajectory

# calc = EAM(potential='FeNiCr_pt.alloy')
calc = EAM(potential='Fe-Ni.eam.alloy')

a = 3.5574

fcc = FaceCenteredCubic('Ni', directions=[[1,0,0],[0,1,0],[0,0,1]], latticeconstant=a)

# slab1 = surface(fcc, (1,1,1), 36)
slab1 = fcc111('Ni', size=(20,20,30), a=a, orthogonal=True)
# slab1 = slab1*(5,5,1)
slab1.set_pbc([1,1,1])

slab2 = fcc111('Ni', size=(20,20,29), a=a, orthogonal=True)
#slab2 = surface(fcc, (1,1,1), 35)
#slab2 = slab2*(5,5,1)
slab2.set_pbc([1,1,1])

print range(421,431)+range(442,451)+range(462,470)+range(483,490)+range(503,509)+range(524,529)+range(544,548)+range(565,567)+[586]
slab3 = slab1.copy()
del slab3[range(421,431)+range(442,451)+range(462,470)+range(483,490)+range(503,509)+range(524,529)+range(544,548)+range(565,567)+[586]]
#slab2 = slab1.copy()
#del slab2[[0,1]]

view(slab1)
write('slab1-Ni.pdb',slab1)
write('slab2-Ni.pdb',slab2)
write('slab3-Ni.pdb',slab3)

#view(slab1)
#view(slab2)
#view(slab3)

slab1.set_calculator(calc)

#opt = BFGS(slab1, trajectory='slab1.traj')
#opt.run(fmax=0.02)

#write('slab1.pdb', slab1)

p1 = slab1.get_potential_energy() / slab1.get_number_of_atoms()
v1 = slab1.get_volume() / slab1.get_number_of_atoms()
a1 = slab1.get_volume()/slab1.get_cell()[2][2]
print p1, v1, a1


#opt = BFGS(slab2, trajectory='slab2.traj')
#opt.run(fmax=0.02)

#write('slab2.pdb', slab2)

slab2.set_calculator(calc)
p2 = slab2.get_potential_energy() / slab2.get_number_of_atoms()
v2 = slab2.get_volume() / slab2.get_number_of_atoms()
a2 = slab2.get_volume()/slab2.get_cell()[2][2]
print p2, v2, a2

print p2-p1, slab1.get_number_of_atoms(), (p2-p1)*96000/6.022e23*slab1.get_number_of_atoms()/a1*1.0e20

slab3.set_calculator(calc)

T = 10
dyn = Langevin(slab3, 5 * units.fs, T * units.kB, 0.002)

traj = Trajectory('md-Ni.traj', 'w', slab3)
dyn.attach(traj.write, interval=50)

def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

printenergy(slab3)

for i in range(200):
    dyn.run(50)
    printenergy(slab3)


