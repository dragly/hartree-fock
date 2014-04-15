#!/usr/bin/python
import subprocess
import os
import os.path
from sys import argv

if len(argv) < 2:
    raise Exception("Missing states_file argument. Usage: run.py <states_file>")

output_dir = os.path.abspath("runs")
if len(argv) < 3:
    project_id = "tmp"
else:
    project_id = argv[-1]
    
output_dir = os.path.join(output_dir, project_id)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
states_file = argv[1]

output_file = os.path.join(output_dir, os.path.split(states_file)[-1])

current_path = os.path.dirname(os.path.realpath(__file__))
build_path = os.path.join("..", "build")
if not os.path.exists(build_path):
    os.makedirs(build_path)
subprocess.call(["qmake", current_path], cwd=build_path)
subprocess.call(["make"], cwd=build_path)

staterunner_path = os.path.join(build_path, "tools", "staterunner")
lib_path = os.path.join("..", "..", "src")

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path
subprocess.call(["mpirun", "-n", "8", "./staterunner", states_file, output_file], cwd=staterunner_path, env=env)
