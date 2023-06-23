# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:56:30 2023

@author: griffin.kowash
"""

import numpy as np
from braid import Braid

# Copied from old Params.set_from_args method
"""
# Command line parameters
arg_str = "s:c:n:d:a:z:r:p:m:f:v:o:"
arg_list = ['plot', 'plotmode =', 'save', 'verbose', 'outdir =', 'fromfile']
options, args = getopt.getopt(sys.argv[1:], arg_str, arg_list)

for opt, arg in options:
    print(opt, arg)
    if opt in ['-s']:
        self.s = float(arg)
    elif opt in ['-c']:
        self.c = int(arg)
    elif opt in ['-n']:
        self.n = int(arg)
    elif opt in ['-d']:
        self.d = float(arg)
    elif opt in ['-a']:
        self.alpha = float(arg)
    elif opt in ['-z']:
        self.z_max = float(arg)
    elif opt in ['-r']:
        self.resolution = int(arg)
    elif opt in ['-p', '--plot']:
        self.plotting = True
    elif opt in ['-m', '--plotmode ']:
        if arg == 'line':
            self.plot_mode = Params.LINES
        elif arg == 'surf':
            self.plot_mode = Params.SURFACES
        else:
            print('Unknown argument for plot mode provided.')
    elif opt in ['-f', '--save']:
        self.saving = True
    elif opt in ['-v', '--verbose']:
        self.verbose = True
    elif opt in ['-o', '--outdir ']:
        self.output_dir = arg
    elif opt in ['--fromfile']:
        self.set_from_file = True
    else:
        print(f'Unknown option {opt} provided.')
"""



#if __name__ == '__main__':
heart = {'range': (0, 2*np.pi), 'func': lambda t: (4*16*np.sin(t)**3, 4*(13*np.cos(t) - 5*np.cos(2*t) - 2*np.cos(3*t) - np.cos(4*t)), 0*np.sin(t))}
spiral = {'range': (0, 50), 'func': lambda t: (0.1 * t * np.cos(t), 0.1 * t * np.sin(t), t)}
square_root = {'range': (0, 10), 'func': lambda t: (t, 5*t**0.5, 0*t)}
parabola = {'range': (-10, 10), 'func': lambda t: (t, 0.02*10*t**2, 0*t)}
quartic = {'range': (-2.5, 2.5), 'func': lambda t: (t, 0.25 * t**4, 0*t)}
worm = {'range': (0, 16), 'func': lambda t: (0*t, 6*np.sin(t), t)}

path_func = heart

#params = Params()
#params.resolution = 300
#params.set_path_from_function(path_func['func'], path_func['range'], equidistant=True)
#params.set_from_spline
#if params.set_from_file:
#    params.set_path_from_file(params.output_dir + '\\path_data.csv', equidistant=True)
#params.verbose = True
#params.plotting = True
#params.saving = True
#params.plot_mode = Params.SURFACES


#braid = Braid(params)
#braid.params.path.plot_xy_nodes(skip=1, both=False)#max(1, int(params.resolution/100)))
#braid.params.path.plot_node_spacing(equidistant=True)
    

braid = Braid()
braid.set_geometry()
#braid.set_path_from_function(path_func['func'], path_func['range'], 600, equidistant=True)
braid.set_linear_path_between((-10,17,0), (3,44,10), 100)
braid.construct(verbose=False)
braid.plot(linewidth=3)
