
from ase.lattice.cubic import *
from ase.calculators.eam import EAM
from ase.visualize import view, write
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt

calc = EAM(potential='FeNiCr_pt.alloy')
calc.plot()

OPTIONS = np.linspace(0.95, 1.05, 11)
volumes = []
energies = []

a = 3.6
bulk = FaceCenteredCubic('Ni', directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], latticeconstant=a)
bulk=bulk*[10,10,10]
a = 36

for opt in OPTIONS:
    bulk.set_cell([a*opt,a*opt,a*opt], scale_atoms=True)

    bulk.set_calculator(calc)
    p = bulk.get_potential_energy()/bulk.get_number_of_atoms()
    v = bulk.get_volume()/bulk.get_number_of_atoms()

    energies.append(p)
    volumes.append(v)

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
eos.plot('eos.png')
plt.show('eos.png')
