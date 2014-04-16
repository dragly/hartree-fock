#!/usr/bin/python
import subprocess
import os
import os.path
import signal
from sys import argv
from argparse import ArgumentParser

try:
    from sumatra.projects import load_project
    project = load_project()
    output_dir = os.path.abspath(project.data_store.root)
except ImportError:
    output_dir = os.path.abspath("tmp")

parser = ArgumentParser()
parser.add_argument("states_file")
parser.add_argument("project_id", nargs='?', default="tmp")
args = parser.parse_args()

current_path = os.path.dirname(os.path.realpath(__file__))

project_id = args.project_id
    
output_dir = os.path.join(output_dir, project_id)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
states_file = argv[1]

output_file = os.path.join(output_dir, os.path.split(states_file)[-1])

build_path = os.path.abspath(os.path.join(current_path, "..", "..", "..", "build"))
project_path = os.path.abspath(os.path.join(current_path, "..", ".."))

print "Building in ", build_path

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
