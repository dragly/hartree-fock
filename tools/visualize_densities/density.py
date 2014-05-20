# Plot the atoms and the bonds ################################################
from pylab import *
from mayavi import mlab
from os.path import join
import h5py
from argparse import ArgumentParser
from glob import glob
import os
import os.path

parser = ArgumentParser()
parser.add_argument("results_path")
parser.add_argument("--id", nargs='?', default="tmp")
args = parser.parse_args()

output_dir = "tmp"

if args.id != "tmp":
    try:
        from sumatra.projects import load_project
        output_dir = os.path.join(os.path.abspath(load_project().data_store.root), args.id)
    except ImportError:
        pass
    
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def draw_atoms(atoms, atom_meta):
    n_electrons = 0
    counter = 0
    for atom in atoms:
        if atom_meta[counter]["type"] == 1:
            color = (1, 1, 1)
        elif atom_meta[counter]["type"] == 8:
            color = (1, 0, 0)
        else:
            color = (1, 1, 0)
        mlab.points3d(atom[0], atom[1], atom[2],
                          scale_factor=0.3,
                          resolution=20,
                          color=color,
                          scale_mode='none')
        counter += 1
    return n_electrons

density_file_name = os.path.join(args.results_path, "density.h5")

atoms_data_file = h5py.File(join(args.results_path, "results.h5"))
atom_meta = atoms_data_file.get("atomMeta")[:]
atoms = atoms_data_file.get("state")[:]
atoms_data_file.close()

density_file = h5py.File(density_file_name)
data = density_file.get("dataset")[:]
density_file.close()

print "Data min,max:",data.min(),data.max()

X,Y,Z = mgrid[-3:3:1j*data.shape[0], -3:3:1j*data.shape[1], -3:3:1j*data.shape[2]]

mlab.figure(2, bgcolor=(0, 0, 0), size=(1280, 720))
mlab.clf()
n_electrons = draw_atoms(atoms, atom_meta)
data_max_min_diff = (data.max() - data.min())
levels = [0.0003, 0.008]
contours = []
for level in levels:
    contours.append(data.min() + level * data_max_min_diff)
contours = [0.005, 0.05, 0.1]
iso = mlab.contour3d(X, Y, Z, data, vmin=contours[0], vmax=contours[-1], opacity=0.5, contours=contours)

mlab.savefig(os.path.join(output_dir, "density.png"))
mlab.savefig(os.path.join(output_dir, "density.x3d"))
