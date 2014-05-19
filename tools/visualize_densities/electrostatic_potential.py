# Plot the atoms and the bonds ################################################
from pylab import *
from mayavi import mlab
from os.path import join
import h5py

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

path_name = "/home/svenni/Dropbox/projects/programming/hartree-fock/build-hartree-fock-Desktop_Qt_5_2_1_GCC_64bit-Release/app"

atoms_data_file = h5py.File(join(path_name, "results.h5"))
atom_meta = atoms_data_file.get("atomMeta")[:]
atoms = atoms_data_file.get("state")[:]
atoms_data_file.close()

file_name = "electrostatic_potential.h5"
density_file = h5py.File(join(path_name, file_name))
data = density_file.get("dataset")[:]
density_file.close()

print "Data min,max:",data.min(),data.max()

X,Y,Z = mgrid[-5:5:1j*data.shape[0], -5:5:1j*data.shape[1], -5:5:1j*data.shape[2]]

mlab.figure(2, bgcolor=(0, 0, 0), size=(1280, 720))
mlab.clf()
n_electrons = draw_atoms(atoms, atom_meta)
data_max_min_diff = (data.max() - data.min())
levels = [0.0003, 0.008]
contours = []
for level in levels:
    contours.append(data.min() + level * data_max_min_diff)
contours = [-0.01, 0.5]
iso = mlab.contour3d(X, Y, Z, data, vmin=contours[0], vmax=contours[-1], opacity=0.5, contours=contours)

mlab.savefig("ch4-volume.x3d")
