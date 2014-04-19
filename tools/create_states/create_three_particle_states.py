import h5py
import numpy
import os
from numpy import dtype, zeros, linspace, pi, cos, sin, sqrt
import yaml
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("config_filename")
parser.add_argument("--id", nargs='?', default="tmp")
args = parser.parse_args()

output_dir = os.path.abspath("tmp")

if args.id != "tmp":
    try:
        from sumatra.projects import load_project
        output_dir = os.path.join(os.path.abspath(load_project().data_store.root), args.id)
    except ImportError:
        pass

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

config_file = open(args.config_filename, "r")
config = yaml.load(config_file)

f = h5py.File("states.h5", "w")

atomType = dtype([("x", float), ("y", float), ("z", float)])
atoms = zeros(3, dtype=atomType)

atomMetaType = dtype([("type", int), ("basisName", "S64")])
atomMeta = zeros(3, dtype=atomMetaType)

atomMeta[0]["type"] = config["type1"]
atomMeta[0]["basisName"] = config["basisName"]
atomMeta[1]["type"] = config["type2"]
atomMeta[1]["basisName"] = config["basisName"]
atomMeta[2]["type"] = config["type3"]
atomMeta[2]["basisName"] = config["basisName"]

dataset2 = f.create_dataset("atomMeta", data=atomMeta)

angles = linspace(config["angleMin"], config["angleMax"], config["angleCount"])
r12s = linspace(config["r12Min"], config["r12Max"], config["r12Count"])
r13s = linspace(config["r13Min"], config["r13Max"], config["r13Count"])

dataset2.attrs["description"] = "H2O with variation of angles and distances"
dataset2.attrs["angleMin"] = angles.min()
dataset2.attrs["angleMax"] = angles.max()
dataset2.attrs["r12Min"] = r12s.min()
dataset2.attrs["r12Max"] = r12s.max()
dataset2.attrs["r13Min"] = r13s.min()
dataset2.attrs["r13Max"] = r13s.max()

statesGroup = f.create_group("states")

skipped = 0

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