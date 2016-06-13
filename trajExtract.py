from ase import Atoms
from ase.visualize import view
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.lattice.surface import fcc111,add_adsorbate
from ase.io import write
from ase.io.trajectory import PickleTrajectory
import io
import os

#Date: 4/17/15

#sets up the slab and molecule
slab = fcc111('Pd', size=(3,3,2), vacuum=10)
molecule = Atoms('CO', [(0., 0., 0.), (0., 0., 1.13)])
add_adsorbate(slab, molecule, 1.7, 'ontop')

#sets up the calculator
calc = EMT()
slab.set_calculator(calc)

#runs optimization
dyn = QuasiNewton(slab, trajectory='CO-Pd-ontop.traj')
dyn.run(fmax=0.05)

#writes the trajectory file in .xyz format
traj = PickleTrajectory("CO-Pd-ontop.traj")
nsteps = dyn.get_number_of_steps()
string = 'structure'

path = os.getcwd()
path = path + '/scratch'
if not os.path.exists(path): os.makedirs(path)

outFileName = 'trajectory.xyz'
for i in range(0, nsteps+1):
	atoms = traj[i]
  	string = 'structure%03d' % (i,) +'.xyz'
	outStruct = os.path.join(path, string)
  	write(outStruct, atoms)
#combines all optimization steps in one trajectory file
	inFile = open (os.path.join(path, 'structure%03d' % (i,) +'.xyz'), 'r')
	fileStr = inFile.read()
	outFile = open(outFileName, 'a')
	outFile.write(fileStr)

#writes the final optimized structure
write('CO-Pd-ontop.xyz', slab)
