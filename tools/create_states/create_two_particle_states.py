import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin
from sys import argv
import os, os.path
import yaml
from argparse import ArgumentParser

try:
    from sumatra.projects import load_project
    project = load_project()
    output_dir = os.path.abspath(project.data_store.root)
except ImportError:
    output_dir = os.path.abspath("tmp")

parser = ArgumentParser()
parser.add_argument("config_filename")
parser.add_argument("project_id", nargs='?', default="tmp")
args = parser.parse_args()

config_file = open(args.config_filename, "r")
output_dir = os.path.join(output_dir, args.project_id)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

config = yaml.load(config_file)

atom_type_1 = config["type1"]
atom_type_2 = config["type2"]
basis_name = config["basisName"]

r12_min = config["r12Min"]
r12_max = config["r12Max"]
configuration_count = config["n"]
r12_ground_state = config["r12GroundState"]

#f = h5py.File("states_atom_" + str(atom_type_1) + "_atom_" + str(atom_type_2) + "_basis_" + basis_name + ".h5", "w")
file_name = "states_atom_%d_atom_%d_basis_%s.h5" % (atom_type_1, atom_type_2, basis_name)
file_path = os.path.join(output_dir, file_name)
f = h5py.File(file_path, "w")

atom_type = dtype([("x", float), ("y", float), ("z", float)])
atoms = zeros(2, dtype=atom_type)

atom_metaType = dtype([("type", int), ("basisName", "S64")])
atom_meta = zeros(2, dtype=atom_metaType)

atom_meta[0]["type"] = atom_type_1
atom_meta[0]["basisName"] = basis_name
atom_meta[1]["type"] = atom_type_2
atom_meta[1]["basisName"] = basis_name

atom_meta_dataset = f.create_dataset("atomMeta", data = atom_meta)
atom_meta_dataset.attrs["r12Min"] = r12_min
atom_meta_dataset.attrs["r12Max"] = r12_max

r12s = linspace(r12_min, r12_max, configuration_count)

if r12_ground_state > 0.0:
    atoms[0]["x"] = -r12_ground_state / 2
    atoms[0]["y"] = 0.0
    atoms[0]["z"] = 0.0
    
    atoms[1]["x"] = r12_ground_state / 2
    atoms[1]["y"] = 0.0
    atoms[1]["z"] = 0.0
    dataset = f.create_dataset("groundState", data=atoms)
    dataset.attrs["r12"] = r12_ground_state

stateCounter = 0
states_group = f.create_group("states")
for j in range(len(r12s)):  
    atoms[0]["x"] = -r12s[j] / 2
    atoms[0]["y"] = 0.0
    atoms[0]["z"] = 0.0
    
    atoms[1]["x"] = r12s[j] / 2
    atoms[1]["y"] = 0.0
    atoms[1]["z"] = 0.0
            
    dataset = states_group.create_dataset("state%010d" % stateCounter, data=atoms)
    dataset.attrs["r12"] = r12s[j]
    stateCounter += 1
    
f.close()
print "States saved to:\n", os.path.abspath(file_path)