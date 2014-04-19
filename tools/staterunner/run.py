#!/usr/bin/python
import subprocess
import os
import os.path
import signal
from sys import argv
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("states_file")
parser.add_argument("--id", nargs='?', default="tmp")
args = parser.parse_args()

output_dir = os.path.abspath("tmp")

if args.id != "tmp":
    try:
        from sumatra.projects import load_project
        output_dir = os.path.join(os.path.abspath(load_project().data_store.root), args.id)
    except ImportError:
        pass

current_path = os.path.dirname(os.path.realpath(__file__))

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
states_file = argv[1]

output_file = os.path.join(output_dir, os.path.split(states_file)[-1])

build_path = os.path.abspath(os.path.join(current_path, "..", "..", "..", "build"))
project_path = os.path.abspath(os.path.join(current_path, "..", ".."))

print "Building in:\n", build_path

if not os.path.exists(build_path):
    os.makedirs(build_path)
subprocess.call(["qmake", project_path, "CONFIG+=nogui"], cwd=build_path)
subprocess.call(["make"], cwd=build_path)

staterunner_path = os.path.join(build_path, "tools", "staterunner")
lib_path = os.path.join("..", "..", "src")

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path
#proc = subprocess.call(["./staterunner", states_file, output_file], cwd=staterunner_path, env=env)
proc = subprocess.call(["mpirun", "-n", "8", "./staterunner", states_file, output_file], cwd=staterunner_path, env=env)

print "Results saved to:\n", output_dir + "/*"