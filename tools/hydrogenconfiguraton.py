import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin
f = h5py.File("hydrogenstates.h5", "w")

atomType = dtype([("x", float), ("y", float), ("z", float)])
atoms = zeros(2, dtype=atomType)

atomMetaType = dtype([("type", int), ("basisName", "S64")])
atomMeta = zeros(2, dtype=atomMetaType)

atomMeta[0]["type"] = 1
atomMeta[0]["basisName"] = "6-311G"
atomMeta[1]["type"] = 1
atomMeta[1]["basisName"] = "6-311G"

f.attrs["atomMetaInformation"] = atomMeta
#f.create_dataset("atomMetaInformation", data=atomMeta)

r12s = linspace(1.0, 10.0, 200)

stateCounter = 0
for j in range(len(r12s)):  
    atoms[0]["x"] = -r12s[j] / 2
    atoms[0]["y"] = 0.0
    atoms[0]["z"] = 0.0
    
    atoms[1]["x"] = r12s[j] / 2
    atoms[1]["y"] = 0.0
    atoms[1]["z"] = 0.0
            
    dataset = f.create_dataset("state%010d" % stateCounter, data=atoms)
    stateCounter += 1
    
f.close()