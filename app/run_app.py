#!/usr/bin/python
import subprocess
import os
import os.path
import signal
from sys import argv
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("config_file")
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
    
config_file = os.path.abspath(argv[1])

output_file = os.path.join(output_dir, os.path.split(config_file)[-1])

build_path = os.path.abspath(os.path.join(current_path, "..", "..", "build-app"))
project_path = os.path.abspath(os.path.join(current_path, ".."))

print "Building in:\n", build_path

if not os.path.exists(build_path):
    os.makedirs(build_path)
subprocess.call(["qmake", project_path, "CONFIG+=nogui CONFIG+=notests CONFIG+=notools"], cwd=build_path)
subprocess.call(["make", "-j", "8"], cwd=build_path)

app_path = os.path.join(build_path, "app")
lib_path = os.path.join(build_path, "src")

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path
#proc = subprocess.call(["./staterunner", states_file, output_file], cwd=staterunner_path, env=env)
run_argument = ["./hartree-fock", config_file]
print " ".join(run_argument)
proc = subprocess.call(run_argument, cwd=app_path, env=env)
