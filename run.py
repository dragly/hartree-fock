#!/usr/bin/python
import subprocess
import os
import os.path
import sys

if len(sys.argv) < 2:
    raise Exception("Missing states_file argument. Usage: run.py <states_file>")    
    
states_file = sys.argv[1]

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
subprocess.call(["./staterunner", states_file], cwd=staterunner_path, env=env)
