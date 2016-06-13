
# Mina Jafari, August 2015
# Purpose: To calculate the forces of final structue of each image in an NEB run. 
# Run this in neb-scratch directory. Define slab based on the structure.

# The vasp path should be set in the shell terminal.

import sys
sys.path.append('/export/zimmerman/paulzim/ase')
from math import sqrt
import os, fnmatch
from ase import Atoms
import io 
from ase.io import read
from subprocess import call
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.lattice.surface import fcc110, add_adsorbate

slab = fcc110('Cu', size=(3,3,4), vacuum=10)
molecule = Atoms('OHH')
add_adsorbate(slab, molecule, 1.7, 'ontop')

path = os.getcwd()

FH = open('finalForces', 'a')
for file in sorted(os.listdir(path)):
    if fnmatch.fnmatch(file, '*.xyz'):
        fName = file
        slabatoms = read(fName)
        slab.set_positions(slabatoms.get_positions())
#    calc = Vasp(xc='PBE',lreal='Auto',kpts=[1,1,1],ismear=1,sigma=0.2,algo='fast',istart=0,npar=8,encut=300)
#        calc = EMT()
        slab.set_calculator(calc)
        f = slab.get_forces()
        fm = sqrt((f**2).sum(axis=1).max()) 
        FH.write(str(file))
        FH.write('\n')
        FH.write(str(fm))
        FH.write('\n')
        FH.write('\n')

FH.close()
