# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:30:29 2023

@author: griffin.kowash
"""

from geomdl import BSpline, convert
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np
import time

# default settings
resolution = 300
plotting = False
spline_file = 'C:\\Users\\griffin.kowash\\AppData\\Local\\Temp\\braids\\braid_data\\spline_data.csv'
output_file = 'C:\\Users\\griffin.kowash\\AppData\\Local\\Temp\\braids\\braid_data\\path_data.csv'

# process command line arguments
if len(sys.argv) > 1:
    arg_str = "d:o:r:"
    arg_list = ['data =', 'output =', 'resolution =', 'plot']
    options, args = getopt.getopt(sys.argv[1:], arg_str, arg_list)
    
    for opt, arg in options:
        if opt in ['-d', '--data']:
            pass
            #spline_file = arg
        elif opt in ['-o', '--output']:
            pass
            #output_file = arg
        elif opt in ['-r', '--resolution']:
            resolution = int(arg)
        elif opt in ['--plot']:
            plotting = True

# read spline data and convert to points
with open(spline_file, 'r') as file:
    lines = file.readlines()

degree = int(lines[0].split(': ')[1][:-1])
knot_vector = [float(knot) for knot in lines[1].split(': ')[1][:-1].split(', ')]

control_points = []
for line in lines[2:]:
    x, y, z = line.split('\n')[0].split(', ')
    control_points.append((1000*float(x), 1000*float(y), 1000*float(z)))

print(degree)
print(knot_vector)
print(control_points)

crv = BSpline.Curve()
crv.degree = degree
crv.ctrlpts = control_points
crv.knotvector = knot_vector

crv_rat = convert.bspline_to_nurbs(crv)

crv_rat.sample_size = resolution
points = np.array(crv_rat.evalpts)


# save file and optionally plot
np.savetxt(output_file, points)
  
if plotting:  
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(points[:, 0], points[:, 1], points[:, 2])
    
    control_points = np.array(control_points)
    ax.scatter(control_points[:, 0], control_points[:, 1], control_points[:, 2])