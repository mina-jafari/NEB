import sys
sys.path.append('/export/zimmerman/mjafari/ase')
sys.path.append('/export/zimmerman/paulzim/vasp.5.3/exe/gollummpi')

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

# Read initial and final states:
initial = read('initial.traj')
final = read('final.traj')
numOfImages = 9

# Set calculators:
images = [initial]
constraint = FixAtoms(mask=[atom.tag > 2 for atom in initial])
calc = Vasp(xc='PBE',lreal='Auto',kpts=[1,1,1],ismear=1,sigma=0.2,algo='fast',istart=0,npar=8,encut=300)
n = size // numOfImages
#j = rank * numOfImages // size 
j = ( ((2 * rank) + 2) // n) - 1
#calc = EMT()
for i in range(0, numOfImages): #determines the number of nodes
    image = initial.copy()
    if i == j:
        image.set_calculator(calc)
    image.set_constraint(constraint)
    images.append(image)
images.append(final)

neb = NEB(images, climb=True, parallel=True)
neb.interpolate()

dyn = BFGS(neb, logfile="logFile") 
#if rank % (size // ) == 0:
#        traj = PickleTrajectory('neb%d.traj' % j, 'w', images[j], master=True)
#        qn.attach(traj)

for i in range(0, numOfImages):
    dyn.attach(PickleTrajectory('neb-%d.traj' % i, 'w', images[i]), master=True)
dyn.run(fmax=0.014)

#writes the coordinates for each image in NEB path in .xyz format
string = 'structure'

path = os.getcwd()
path = path + '/neb-scratch'
if not os.path.exists(path): os.makedirs(path)

outFileName = 'NEB-trajectory.xyz'
if os.path.exists(outFileName):
    os.remove(outFileName)

for i in range(0, neb.nimages):
    string = 'structure%03d' % (i,) +'.xyz'
    outStruct = os.path.join(path, string)
    write(outStruct, images[i])
    f = open(outStruct, 'a')
    f.write(str(images[i].get_potential_energy() - images[0].get_potential_energy()))
    f.close()

for i in range(0, neb.nimages):
#combines all images into one file with energies written for each structure
    inFile = open(os.path.join(path, 'structure%03d' % (i,) +'.xyz'), 'r')
    linesInFile = list()
    for line in inFile:
        linesInFile.append(line)
    numOfLines = len(linesInFile)
    outFile = open(outFileName, 'a')
    outFile.write(linesInFile[0])
    outFile.write(linesInFile[numOfLines-1])
    for i in range(1, numOfLines-1):
        outFile.write(linesInFile[i])

#writes the trajectory file in .xyz format
if os.path.exists("gradient"):
    os.remove("gradient")
nsteps = dyn.get_number_of_steps()
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


fh = open("gradient", 'a')
fh.write('\n\n\n')
fh.write("Total gradients = ")
fh.write(str((nsteps+1) * (neb.nimages-2)))
fh.close()
