import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin
f = h5py.File("hydrogenstates.h5", "r")

atomsMeta = f.attrs["atomMetaInformation"]

out_file = open("hydrogenstates.xyz", "w")

for stateName in f:
    if "state" in stateName:
        state = f.get(stateName)
        # Write header
        out_file.write(str(len(state)) + "\n")
        out_file.write("\n")
        for i in range(len(atomsMeta)):
            out_file.write("%d %.10e %.10e %.10e\n" % (atomsMeta[i]["type"], state[i]["x"], state[i]["y"], state[i]["z"]))
            
out_file.close()