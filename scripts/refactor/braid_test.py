# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:56:30 2023

@author: griffin.kowash
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from braid import Braid
from spline import Spline

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

path_func = parabola

#spline_path = r'C:\Users\griffin.kowash\AppData\Local\Temp\braids\braid_data\spline.json'
#with open(spline_path, 'r') as spline_file:
#    spline_dict = json.load(spline_file)

"""
### Vance Fig. 8 ###
geo = {
    's': 1,
    'c': 42,
    'n': 10,
    'd': 0.015748,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0,#0.06,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 0,#5,
    'tanh_c': 0#0.5
    }
"""

### Yazaki ###
geo = {
    's': 10.08,
    'c': 42,
    'n': 9,
    'd': 0.16,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.08,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 5,
    'tanh_c': 0.5
    }


braid = Braid()
braid.set_geometry(geo)
#braid.set_path_from_spline(spline_dict, 600, equidistant=False)
#braid.set_path_from_function(path_func['func'], path_func['range'], 200)
#braid.set_linear_path_between((0,0,-34), (0,0,112.04), 2300)
braid.set_linear_path_between((0,0,0), (0,0,60), 200)
#braid.set_linear_path_between((0,0,-34), (0,0,136), 500)
braid.construct(mode='surface', verbose=False, testing=False)
braid.plot(linewidth=0.5, mode='lines')
braid.save(r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\braids\Vance_k80_surface', fmt='default', mode='surface')



def calc_params(geo):
    W = 4 * np.pi * geo['s'] / geo['c'] * np.cos((90 - geo['alpha_deg']) * np.pi / 180)
    F = (geo['n']) * geo['d'] / W  # note that geo['n'] is currently being treated as number of strands + 1
    K = 2*F - F**2
    print(f'Shield coverage (K): {round(K * 100, 1)}%')
    
calc_params(geo)

