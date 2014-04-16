import h5py
from pylab import *
from glob import glob
from sys import argv
import os
import os.path
import time
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("states_file", nargs="?")
parser.add_argument("--parentid")
parser.add_argument("--parentproject")
#parser.add_argument("project_id", nargs='?', default="tmp")
args = parser.parse_args()

output_dir = os.path.abspath("tmp")
project_id = "tmp"

record_run = True

if record_run:
    from sumatra.projects import load_project
    project = load_project()
    output_dir = os.path.abspath(load_project().data_store.root)
    record = project.new_record(main_file=__file__)
    project_id = record.label
    start_time = time.time()

output_dir = os.path.join(output_dir, project_id)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if args.parentid and args.parentproject:
    from sumatra.projects import load_project
    project = load_project()
    parent_record = project.record_store.get(args.parentproject, args.parentid)
    parent_output_data = parent_record.output_data
    states_files = []
    for key in parent_output_data:
        states_files.append(os.path.join(parent_record.datastore.root, key.path))
        record.input_data.append(key)
else:    
    states_file = args.states_file
    if os.path.isdir(states_file):
        states_files = glob(states_file + "/*.h5")
    else:
        states_files = glob(states_file)

output_file = os.path.join(output_dir, "two_particle_plot.pdf")

energy_min = inf
energy_max = -inf

energies = []
r12s = []

for statesFile in states_files:
    f = h5py.File(statesFile, "r")
    atomsMeta = f.get("atomMeta")
    if len(atomsMeta) != 2:
        raise Exception("Wrong number of atoms in atomsMeta. Found %d, should be 2." % len(atomsMeta))
    plot_title = "Distance between atom %d and %d" % (atomsMeta[0]["type"], atomsMeta[1]["type"])
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

#figure()
title(plot_title)
plot(r12s, energies)
xlabel(r"$r$")
ylabel(r"$E$")
savefig(output_file)
savefig(output_file + ".png")

print "Results saved to", output_file

if record_run:
    record.duration = time.time() - start_time
    record.output_data = record.datastore.find_new_data(record.timestamp)
    project.add_record(record)
    project.save()