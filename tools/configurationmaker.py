import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin
f = h5py.File("states.h5", "w")

atomType = dtype([("x", float), ("y", float), ("z", float)])
atoms = zeros(3, dtype=atomType)

atomMetaType = dtype([("type", int), ("basisName", "S64")])
atomMeta = zeros(3, dtype=atomMetaType)

atomMeta[0]["type"] = 8
atomMeta[0]["basisName"] = "3-21G"
atomMeta[1]["type"] = 1
atomMeta[1]["basisName"] = "3-21G"
atomMeta[2]["type"] = 1
atomMeta[2]["basisName"] = "3-21G"

dataset2 = f.create_dataset("atomMeta", data=atomMeta)

angles = linspace(pi/3, pi, 10)
r12s = linspace(1.0, 6.0, 10)
r13s = linspace(1.0, 6.0, 10)

dataset2.attrs["description"] = "H2O with variation of angles and distances"
dataset2.attrs["angleMin"] = angles.min()
dataset2.attrs["angleMax"] = angles.max()
dataset2.attrs["r12Min"] = r12s.min()
dataset2.attrs["r12Max"] = r12s.max()
dataset2.attrs["r13Min"] = r13s.min()
dataset2.attrs["r13Max"] = r13s.max()

statesGroup = f.create_group("states")

stateCounter = 0
for j in range(len(r12s)):  
    for k in range(len(r13s)):
        for i in range(len(angles)):
            atoms[0]["x"] = 0.0
            atoms[0]["y"] = 0.0
            atoms[0]["z"] = 0.0
            
            atoms[1]["x"] = r12s[j]
            atoms[1]["y"] = 0.0
            atoms[1]["z"] = 0.0
            
            atoms[2]["x"] = cos(angles[i]) * r13s[k]
            atoms[2]["y"] = sin(angles[i]) * r13s[k]
            atoms[2]["z"] = 0.0
            
            dataset = statesGroup.create_dataset("state%010d" % stateCounter, data=atoms)

            dataset.attrs["angle"] = angles[i]
            dataset.attrs["r12"] = r12s[j]
            dataset.attrs["r13"] = r13s[k]
            
            stateCounter += 1
    
f.close()