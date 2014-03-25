# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 17:35:33 2014

@author: svenni
"""

from pylab import imshow, figure, loadtxt, show

def rescale(value, value_min, value_max):
    return (value - value_min) / (value_max - value_min) * 0.8 + 0.1

def rescale_inverse(value, value_min, value_max):
    return (value - 0.1) / 0.8 * (value_max - value_min) + value_min

data = loadtxt("/home/svenni/Dropbox/projects/programming/build-emdee-Desktop_Qt_5_2_0_with_GDB-Release/tools/fann-plot-test/testdata5.444")

energy_min=-75.5556102338
energy_max=-74.1449902271

data = rescale_inverse(data, energy_min, energy_max)

figure()
imshow(data, interpolation="nearest", origin="lower", vmin=energy_min, vmax=energy_max)
show()
