import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin, arctan, arctan2, sqrt, meshgrid, array
from pylab import imshow, plot, figure, subplot, title
from glob import glob
from scipy.interpolate import griddata
from collections import OrderedDict

plots = []

for statesFile in glob("/home/svenni/Dropbox/studies/master/code/hartree-fock/build-hartree-fock-Desktop_Qt_5_2_0_with_GDB-Release/tools/staterunner/results.h5.*"):
    #f = h5py.File("/home/svenni/Dropbox/studies/master/code/hartree-fock/build-hartree-fock-Desktop_Qt_5_2_0_with_GDB-Release/tools/staterunner/results.h5.0000", "r")
    f = h5py.File(statesFile, "r")
    atomsMeta = f.attrs["atomMetaInformation"]
    print statesFile
    for stateName in f:
        atoms = f.get(stateName)
        r12 = atoms[1]["x"] - atoms[0]["x"]
        plots.append([r12, atoms.attrs["energy"]])
    f.close()

plots = array(plots)

fig = figure()
plot(plots[:,0], plots[:,1], 'o')
show()