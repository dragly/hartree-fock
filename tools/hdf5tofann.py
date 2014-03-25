import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin, arctan, arctan2, sqrt, meshgrid, array, inf
from pylab import imshow, plot, figure, subplot, title
from glob import glob
from scipy.interpolate import griddata
from collections import OrderedDict

def rescale(value, value_min, value_max):
    return (value - value_min) / (value_max - value_min) * 0.8 + 0.1

file_name = "/home/svenni/Dropbox/studies/master/code/hartree-fock/build-hartree-fock-Desktop_Qt_5_2_0_with_GDB-Release/tools/staterunner/results.h5.*"

train_file = open("train.fann", "w")
test_file = open("test.fann", "w")

n_states_total = 0
n_parameters = 3
n_outputs = 1

energy_min = inf
energy_max = -inf
r12_min = inf
r12_max = -inf
r13_min = inf
r13_max = -inf
angle_min = inf
angle_max = -inf

for statesFile in glob(file_name):
    f = h5py.File(statesFile, "r")
    atomsMeta = f.get("atomMeta")
    r12_min = atomsMeta.attrs["r12Min"]
    r12_max = atomsMeta.attrs["r12Max"]
    r13_min = atomsMeta.attrs["r13Min"]
    r13_max = atomsMeta.attrs["r13Max"]
    angle_min = atomsMeta.attrs["angleMin"]
    angle_max = atomsMeta.attrs["angleMax"]
    states = f.get("/states")
    n_states_total += len(states)
    for stateName in states:
        atoms = states.get(stateName)
        energy_min = min(energy_min, atoms.attrs["energy"])
        energy_max = max(energy_max, atoms.attrs["energy"])
        
    f.close()

n_states_train = int(n_states_total * 0.8)
n_states_test = n_states_total - n_states_train

print n_states_train, n_states_test
    
train_file.write("%d %d %d\n\n" % (n_states_train, n_parameters, n_outputs))
test_file.write("%d %d %d\n\n" % (n_states_test, n_parameters, n_outputs))

state_counter = 0

for statesFile in glob(file_name):
    f = h5py.File(statesFile, "r")
    atomsMeta = f.get("atomMeta")
    states = f.get("/states")
    for stateName in states:
        if state_counter < n_states_train:
            target_file = train_file
        else:
            target_file = test_file
            
        atoms = states.get(stateName)
        r12 = rescale(atoms.attrs["r12"], r12_min, r12_max)
        r13 = rescale(atoms.attrs["r13"], r13_min, r13_max)
        angle = rescale(atoms.attrs["angle"], angle_min, angle_max)
        energy = rescale(atoms.attrs["energy"], energy_min, energy_max)
        
        target_file.write("%.10f %.10f %.10f\n\n" % (r12, r13, angle))
        target_file.write("%.10f\n\n" % (energy))
        state_counter += 1
        
    f.close()
    
train_file.close()
test_file.close()