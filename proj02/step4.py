from ase.lattice.cubic import *
from ase.lattice.surface import surface
from ase.calculators.eam import EAM
from ase.visualize import view, write
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
import elastic
from ase.optimize import BFGS
from ase.build import fcc111
from ase.constraints import FixAtoms

from ase.md.langevin import Langevin
from ase import units
from ase.io.trajectory import Trajectory

calc = EAM(potential='Fe-Ni.eam.alloy')

a = 3.5574

slab1 = fcc111('Ni', size=(20,20,30), a=a, orthogonal=True)
slab1.set_pbc([1,1,1])

slab2 = fcc111('Ni', size=(20,20,29), a=a, orthogonal=True)
slab2.set_pbc([1,1,1])

print range(421,431)+range(442,451)+range(462,470)+range(483,490)+range(503,509)+range(524,529)+range(544,548)+range(565,568)+range(585,587)+[606]
slab3 = slab1.copy()
del slab3[range(5421,5436)]

# del slab3[range(1,15)]

mask = [ atom.position[0] < 1.5 or atom.position[1] < 1.0 or atom.position[2] < 0.1 for atom in slab3]
constraint = FixAtoms(mask=mask)
slab3.set_constraint(constraint)

view(slab3)
write('slab3-Ni-disl.pdb',slab3)

slab3.set_calculator(calc)

#opt = BFGS(slab3, trajectory='slab3.traj')
# opt.run(fmax=0.02)

n = slab3.get_number_of_atoms()
p = slab3.get_potential_energy() / n
v = slab3.get_volume() / n

print n, p, v

T = 10
dyn = Langevin(slab3, 5 * units.fs, T * units.kB, 0.002)

traj = Trajectory('md-Ni-disl.traj', 'w', slab3)
dyn.attach(traj.write, interval=10)

def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

printenergy(slab3)

for i in range(200):
    dyn.run(10)
    printenergy(slab3)