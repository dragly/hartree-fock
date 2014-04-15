import h5py
import numpy
from numpy import dtype, zeros, linspace, pi, cos, sin
from sys import argv
import os, os.path
import yaml

print "Input was", argv

output_dir = "../../runs"

if len(argv) < 2:
    project_id = "unnamed"
else:
    project_id = argv[-1]

output_dir = os.path.join(output_dir, project_id)
    
f = open(os.path.join(output_dir, "test.txt"), "w")
f.write("This is a test that opened " + argv[1])
f.close()
