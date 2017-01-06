from ase.lattice.cubic import *
from ase.lattice.surface import surface
from ase.calculators.eam import EAM
from ase.visualize import view, write
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
import elastic
from ase.optimize import BFGS
from ase.build import fcc111
from ase.lattice import bulk

calc = EAM(potential='FeNiCr_pt.alloy')

a = 3.5574

fcc = FaceCenteredCubic('Fe', directions=[[1,0,0],[0,1,0],[0,0,1]], latticeconstant=a, pbc=True)

fcc.set_calculator(calc)

n = fcc.get_number_of_atoms()
p = fcc.get_potential_energy() / n
v = fcc.get_volume() / n

print n, p, v

slab1 = surface(fcc, (1,1,1), 3)
slab1 = slab1*(1,1,1)

slab1.set_calculator(calc)

n = slab1.get_number_of_atoms()
p = slab1.get_potential_energy() / n
v = slab1.get_volume() / n

print n, p, v

slab2 = fcc111('Fe', size=(4,4,20), a=a, orthogonal=False)
slab2.set_pbc([True, True, True])

view(slab2)

slab2.set_calculator(calc)

n = slab2.get_number_of_atoms()
p = slab2.get_potential_energy() / n
v = slab2.get_volume() / n

print n, p, v

a = 3.5574
a0 = a / np.sqrt(2)
c0 = np.sqrt(8 / 3.0) * a0
# c0 = 1.63*a0
print np.sqrt(8 / 3.0)
hcp = bulk('Fe', 'hcp', a=a0, c=c0)

hcp.set_calculator(calc)

opt = BFGS(hcp, trajectory='hcp.traj')
opt.run(fmax=0.02)

n = hcp.get_number_of_atoms()
p = hcp.get_potential_energy() / n
v = hcp.get_volume() / n

print n, p, v
