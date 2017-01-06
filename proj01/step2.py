from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase import units
from ase.calculators.eam import EAM
from ase.visualize import view, write
from ase.io.trajectory import Trajectory

calc = EAM(potential='FeNiCr_pt.alloy')

size = 5

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol='Fe',
                          size=(size, size, size),
                          pbc=True,
                          latticeconstant=3.5646)

view(atoms)

atoms.set_calculator(calc)

# MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# dyn = VelocityVerlet(atoms, 5 * units.fs)
T = 10
dyn = Langevin(atoms, 5 * units.fs, T * units.kB, 0.002)

traj = Trajectory('md.traj', 'w', atoms)
dyn.attach(traj.write, interval=50)

def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

printenergy(atoms)

for i in range(200):
    dyn.run(50)
    printenergy(atoms)