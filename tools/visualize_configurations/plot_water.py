import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin, arctan, arctan2, sqrt, meshgrid, array, inf
from pylab import imshow, plot, figure, subplot, title
from glob import glob
from scipy.interpolate import griddata
from collections import OrderedDict

plots = {}

energy_min = inf
energy_max = -inf



for statesFile in glob("/home/svenni/Dropbox/studies/master/code/hartree-fock/build-hartree-fock-Desktop_Qt_5_2_1_GCC_64bit-Release/tools/staterunner/results.h5.*"):
    f = h5py.File(statesFile, "r")
    atomsMeta = f.get("atomMeta")
    states = f.get("/states")
    for stateName in states:
        atoms = states.get(stateName)
        #r12 = atoms[1]["x"] - atoms[0]["x"]
        #r13 = sqrt((atoms[2]["x"] - atoms[0]["x"])**2 + (atoms[2]["y"] - atoms[0]["y"])**2)
        #angle = arctan2(atoms[2]["y"], atoms[2]["x"])
        r12 = atoms.attrs["r12"]
        r12_name = "%.4f" % r12
        r13 = atoms.attrs["r13"]
        angle = atoms.attrs["angle"]
        energy = atoms.attrs["energy"]
        if not plots.has_key(r12_name):
            plots[r12_name] = {"angles_r13s": [], "energies": []}
        plots[r12_name]["angles_r13s"].append([angle, r13])
        plots[r12_name]["energies"].append(energy)
        
        r13_min = atomsMeta.attrs["r13Min"]
        r13_max = atomsMeta.attrs["r13Max"]
        
        angle_min = atomsMeta.attrs["angleMin"]
        angle_max = atomsMeta.attrs["angleMax"]
        
        energy_min = min(energy_min, energy)
        energy_max = max(energy_max, energy)
    f.close()

print "Angle min:", angle_min, "max:", angle_max
print "r13 min:", r13_min, "max:", r13_max

angles = linspace(angle_min, angle_max, 100)

r13s = linspace(r13_min, r13_max, 100)
grid_angles, grid_r13s = meshgrid(angles, r13s)
n_plots = len(plots)
n_plots_per_dim = int(sqrt(n_plots) + 1)
plot_counter = 1
plots = OrderedDict(sorted(plots.items()))

print "vmin,vmax: ", energy_min, energy_max

fig = figure(figsize=(20,20))
for r12_name in plots:
    
    ax = fig.add_subplot(n_plots_per_dim, n_plots_per_dim, plot_counter)
    title(r12_name)
    
    values = plots[r12_name]
    
    grid_energies = griddata(array(values["angles_r13s"]), array(values["energies"]), (grid_angles, grid_r13s), method="nearest")
    
    img = imshow(grid_energies, origin="lower", vmin=energy_min, vmax=energy_max,
                 extent=[angle_min,angle_max, r13_min, r13_max], 
                         interpolation="nearest", 
                         aspect=(angle_max - angle_min) / (r13_max - r13_min))
    colorbar()
    
    plot_counter += 1
    