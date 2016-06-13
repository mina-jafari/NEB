import sys
sys.path.append('/export/zimmerman/mjafari/ase')

from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.io import write, read
from ase.neb import NEB
from ase.io.trajectory import PickleTrajectory
import io, os
from ase.parallel import rank, size

numOfImages = 9
f = open('logFile', 'r')
for line in f:
    pass
last = line
split = last.split()
nsteps = int(split[1]) - 1

#writes the trajectory file in .xyz format
if os.path.exists("gradient"):
    os.remove("gradient")

for j in range(0, numOfImages):
    traj = PickleTrajectory('neb-%01d.traj' % (j,))
    path1 = os.getcwd()
    path1 = path1 + '/scratch%01d' % (j,) 
    if not os.path.exists(path1): os.makedirs(path1)
    outFileName1 = 'trajectory%03d' % (j,) + '.xyz'

    fh = open("gradient", 'a')
    fh.write("\n")
    fh.write(str('********neb-%01d.traj' % (j,)))
    fh.write("      Printing information for image %01d **********" % ((j+1),))
    for i in range(0, nsteps+1):
        atoms = traj[i]
        string1 = 'structure%03d' % (i,) +'.xyz'
        outStruct1 = os.path.join(path1, string1)
        if os.path.exists(outStruct1):
            os.remove(outStruct1)
        write(outStruct1, atoms)
        f = open(outStruct1, 'a')
        f.write(str(atoms.get_potential_energy()))
        f.close()

        fh.write('\n')
        fh.write("optimization step %03d" % (i,))
        fh.write("          Energy: ")
        fh.write(str(atoms.get_potential_energy()))
        fh.write("          ")
        fh.write("dE: ")
        fh.write(str(traj[i].get_potential_energy()-traj[i-1].get_potential_energy()))
    fh.write("\n")
    fh.close()

outFileName = 'NEB-traj.xyz'
path = os.getcwd()
for i in range(0, 9):
    dir = path + '/scratch%01d' % (i,)
    os.chdir(dir)
    inFile = open( 'structure%03d' % (nsteps,) +'.xyz', 'r')
    linesInFile = list()
    for line in inFile:
        linesInFile.append(line)
    numOfLines = len(linesInFile)
    os.chdir(path)
    outFile = open(outFileName, 'a')
    outFile.write(linesInFile[0])
    outFile.write(linesInFile[numOfLines-1])
    for j in range(1, numOfLines-1):
        outFile.write(linesInFile[j])
    
fh = open("gradient", 'a')
fh.write('\n\n\n')
fh.write("Total gradients = ")
fh.write(str((nsteps+1) * (11-2)))
fh.close()
