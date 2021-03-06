import h5py
from pylab import *
from glob import glob
from sys import argv
import os
import os.path
import time
from argparse import ArgumentParser
from scipy.optimize import curve_fit
   
def shifted_lennard_jones(x, epsilon, sigma, a):
    return 4 * epsilon * ((sigma / (a + x))**12 - (sigma/(a + x))**6)

parser = ArgumentParser()
parser.add_argument("states_files", nargs="+")
parser.add_argument("--id", default="tmp")
args = parser.parse_args()

output_dir = os.path.abspath("tmp")

if args.id != "tmp":
    try:
        from sumatra.projects import load_project
        output_dir = os.path.join(os.path.abspath(load_project().data_store.root), args.id)
    except ImportError:
        pass

states_files = args.states_files
if len(states_files) == 1:
    states_files = glob(states_files[0])

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_file = os.path.join(output_dir, "two_particle_plot.pdf")

energy_min = inf
energy_max = -inf

energies = []
r12s = []

for statesFile in states_files:
    f = h5py.File(statesFile, "r")
    atomsMeta = f.get("atomMeta")
    energyOffset = atomsMeta.attrs["energyOffset"]
    if len(atomsMeta) != 2:
        raise Exception("Wrong number of atoms in atomsMeta. Found %d, should be 2." % len(atomsMeta))
    states = f.get("/states")
    for stateName in states:
        atoms = states.get(stateName)
        #r12 = atoms[1]["x"] - atoms[0]["x"]
        #r13 = sqrt((atoms[2]["x"] - atoms[0]["x"])**2 + (atoms[2]["y"] - atoms[0]["y"])**2)
        #angle = arctan2(atoms[2]["y"], atoms[2]["x"])
        r12 = atoms.attrs["r12"]
        energy = atoms.attrs["energy"]
        
        r12s.append(r12)
        energies.append(energy)
        
        energy_min = min(energy_min, energy)
        energy_max = max(energy_max, energy)
    f.close()

# Sort r12 and energies by r12
r12s, energies = [list(x) for x in zip(*sorted(zip(r12s, energies), key=lambda pair: pair[0]))]

r12s = array(r12s)
energies = array(energies)

diffs = abs(diff(energies) / diff(r12s))

print "Plotting", len(r12s), "data points."

grid()
#plot(r12s, energies - energyOffset)
plot(r12s, energies)
#popt,errs = curve_fit(shifted_lennard_jones, r12s, energies - energyOffset, p0=(0.01938786,  1.93348366, -0.65246251))
#plot(r12s, shifted_lennard_jones(r12s, *popt))
#plot(r12s[:-1], diffs)
#plot(r12s, 0.13*((1.41/r12s)**12 - 2*(1.41/r12s)**6) - 1.0)
xlabel(r"$r$")
ylabel(r"$E$")
savefig(output_file)
savefig(output_file + ".png")

print "Results saved to", output_file