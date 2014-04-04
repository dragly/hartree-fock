# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 17:35:33 2014

@author: svenni
"""

from pylab import imshow, figure, loadtxt, show, sqrt, colorbar, inf, title
from glob import glob
from os.path import split

def rescale(value, value_min, value_max):
    return (value - value_min) / (value_max - value_min) * 0.8 + 0.1

def rescale_inverse(value, value_min, value_max):
    return (value - 0.1) / 0.8 * (value_max - value_min) + value_min

files = glob("/home/svenni/Dropbox/projects/programming/fann-md/build-fann-md-jimbo-Desktop_Qt_5_2_1_GCC_64bit-Release/fann-plot-test/energies_*.dat")
files = sorted(files)
plot_counter = 1
n_plots_per_dim = int(sqrt(len(files)) + 1)
fig = figure()

value_min = inf
value_max = -inf

for data_file in files:
    data = abs(loadtxt(data_file))
    value_min = min(value_min, data.min())
    value_max = max(value_max, data.max())

for data_file in files:
    ax = fig.add_subplot(n_plots_per_dim, n_plots_per_dim, plot_counter)
    data = abs(loadtxt(data_file))
    
    title(split(data_file)[-1])
    
    imshow(data, interpolation="nearest", origin="lower", vmin=value_min, vmax=value_max)
    colorbar()
    show()
    plot_counter += 1