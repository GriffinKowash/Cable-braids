# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:56:30 2023

@author: griffin.kowash
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from braid_v2 import Braid
from spline import Spline
from scipy.special import ellipk
from scipy.special import ellipe

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


# Same dimensions as Vance Fig. 8; C and N set to 80% coverage
Vance_k80 = {
    's': 10.08,
    'c': 42,
    'n': 9,
    'd': 0.16,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.00,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 12,
    'tanh_c': 0.3
    }

# Test for Simutech Steve's braid
Steve = {
    's': 10.047,
    'c': 48,
    'n': 18,
    'd': 0.127,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.00,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 12,
    'tanh_c': 0.3
    }

# Original validation braid from Yazaki
Yazaki = {
    's': 0.94e-3,
    'ID': 1.64e-3,
    'c': 16,
    'n': 5,
    'd': 0.12e-3,
    'alpha_deg': 19.15,
    'twist_deg': 0,#-7.53,
    'h': 0.06e-3,
    'tanh_m': 6,
    'tanh_c': 0.4,
    'weave_offset': 0,
    'gap_shift_frac': 0
    }

Vance_k85 = {
    's': 10.00 + 0.08,
    'c': 42,
    'n': 10,
    'd': 0.16,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.00,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 12,
    'tanh_c': 0.3
    }

Vance_k95 = {
    's': 10.00 + 0.08,  #mean radius--note that this is typically given as inner radius in literature (and with the variable name "a")
    'c': 52,
    'n': 10,
    'd': 0.16,
    'alpha_deg': 90 - 30,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.16,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 12,
    'tanh_c': 0.3
    }

Vance_k69 = {
    's': 8.00 + 0.16,  #mean radius--note that this is typically given as inner radius in literature (and with the variable name "a")
    'c': 24,
    'n': 11,
    'd': 0.16,
    'alpha_deg': 90 - 20,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.16/2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 12,
    'tanh_c': 0.3
    }

# Braid n.8 from Schippers et al.
Schippers_8 = {
    's': 3.160,  #mean radius
    'c': 24,
    'n': 7 + 1,
    'd': 0.160,
    'alpha_deg': 90 - 38.6,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.080 * 1.0,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 7,
    'tanh_c': 0.5,
    'wire_spread': 1.0,
    'carrier_spread': 1.0
    }

# Braid n.12 from Schippers et al.
Schippers_12 = {
    's': 4.202,  #mean radius
    'c': 32,
    'n': 5 + 1,
    'd': 0.202,
    'alpha_deg': 90 - 35.2,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.202,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.7
    }

# Schippers_12 with radius modified so hole and braid inductance are equal.
Equiductance = {
    's': 4.80 + 0.202,  #mean radius
    'c': 32,
    'n': 6,  #N=5 for physical braid; using N=6 to reach edges of carriers
    'd': 0.202,
    'alpha_deg': 90 - 35.2,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.202,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.6
    }

# Sample 1 from Otin & Isanta 2013
Otin13_1 = {
    's': 0.750 + 0.100,  #mean radius
    'c': 16,
    'n': 4,  #N=3 for physical braid; using N=4 to reach edges of carriers
    'd': 0.100,
    'alpha_deg': 90 - 21.44,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.050,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5
    }

# Sample 2 from Otin & Isanta 2013
Otin13_2 = {
    's': 1.475 + 0.101,  #mean radius
    'c': 16,
    'n': 7,  #N=6 for physical braid; using N=7 to reach edges of carriers
    'd': 0.101,
    'alpha_deg': 90 - 30.08,
    'twist_deg': 0,#-7.53,
    'wave_width': 0.101 / 2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5
    }

# Sample 2 from Otin et al. 2015
Otin15_2 = {
    's': 3.00 + 0.160,  #mean radius
    'c': 24,
    'n': 7,
    'd': 0.160,
    'h': 0.5 * 0.160,
    'alpha_deg': 90 - 32.20,
    'twist_deg': 0,#-7.53,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5,
    'wire_spread': 1.05
    }

# Sample 5 from Otin et al. 2015
Otin15_5 = {
    's': 0.75 + 0.100,  #mean radius
    'c': 16,
    'n': 3,
    'd': 0.100,
    'h': 0.047,
    'alpha_deg': 90 - 21.44,
    'twist_deg': 0,#-7.53,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5,
    'wire_spread': 1.10
    }

# Glenair braid 100-002A1000
Glenair_100002A1000 = {
    's': 12.7 + 0.127,  #mean radius
    'c': 64,
    'n': 12 + 1,
    'd': 0.127,
    'alpha_deg': 90 - 38.8,  #K90=26.6, K95=38.8
    'twist_deg': 0,#-7.53,
    'wave_width': 0.127 / 2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5,
    'skip_factor': 3
    }

# Glenair braid 100-002A125
Glenair_100002A125 = {
    's': 1.6 + 0.127,  #mean radius
    'ID': 3.2,
    'c': 24,
    'n': 5,
    'd': 0.127,
    'h': 0.7 * 0.127,
    'alpha_deg': 12.5,  #K95=12.5 (K90 not physically viable. Some ambiguity due to small size of braid.)
    'twist_deg': 0,#-7.53,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5,
    'wire_spread': 1.10
    }

# Glenair braid 100-002A250
Glenair_100002A250 = {
    's': 3.2 + 0.127,  #mean radius
    'c': 24,
    'n': 16 + 1,
    'd': 0.127,
    'alpha_deg': 90 - 17.6,   #k95=17.6, k90=27.6
    'twist_deg': 0,#-7.53,
    'wave_width': 0.127 / 2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5
    }

# Glenair braid 100-002A500
Glenair_100002A500 = {
    's': 6.35 + 0.127,  #mean radius
    'c': 48,
    'n': 11 + 1,
    'd': 0.127,
    'alpha_deg': 90 - 16,  
    'twist_deg': 0,#-7.53,
    'wave_width': 0.127 / 2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5
    }

# Glenair braid 100-002A750
Glenair_100002A750 = {
    's': 9.55 + 0.127,  #mean radius
    'c': 48,
    'n': 18 + 1,
    'd': 0.127,
    'alpha_deg': 90 - 10,  
    'twist_deg': 0,#-7.53,
    'wave_width': 0.127 / 2,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 1,
    'tanh_c': 0.5,
    'skip_factor': 4
    }

# testing sandbox
test_braid = {
    's': 1,  #mean radius
    'c': 16,
    'n': 9,
    'd': 0.100,
    'alpha_deg': 90 - 20,
    'twist_deg': 0,#-7.53,
    'wave_width': 3,
    'wave_sharpness': 0,#1/3,
    'gap_shift_frac': 0,
    'tanh_m': 7,
    'tanh_c': 0.5,
    'wire_spread': 1.00,
    'carrier_spread': 1.00
    }


geo = Yazaki

braid = Braid()
braid.set_geometry(geo)
#braid.set_path_from_spline(spline_dict, 600, equidistant=False)
#braid.set_path_from_function(path_func['func'], path_func['range'], 200)
#braid.set_linear_path_between((0,0,-34), (0,0,112.04), 2300)
#braid.set_linear_path_between((0,0,-35), (0,0,2), 300)
#braid.set_linear_path_between((0,0,0), (0,0,1), 300)
braid.set_linear_path_between((0,0,-34e-3), (0,0,1.12e-3), 300)
#braid.set_linear_path_between((0,0,-60), (0,0,1.12), 600)
#braid.set_linear_path_between((0,0,-90), (0,0,1.12), 600)
#braid.set_linear_path_between((0,0,-136), (0,0,1.12), 1200)
#braid.set_linear_path_between((0,0,0), (0,0,60), 300)
#braid.set_linear_path_between((0,0,-34), (0,0,136), 500)
braid.construct(verbose=False, testing=False)
braid.plot(linewidth=0.5, mode='lines')
#braid.save(r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\braids\Glenair_100-002A125_body_h70', fmt='default', mode='line')



def calc_params(geo, old=False):
    if old:
        alpha = 90 - geo['alpha_deg']
    else:
        alpha = geo['alpha_deg']
    
    W = 4 * np.pi * (geo['s']) / geo['c'] * np.cos(alpha * np.pi / 180)
    F = (geo['n']) * geo['d'] / W
    K = 2*F - F**2
    
    u = 1.257E-6
    Dm = 2 * geo['s']  # this will break if s is set to inner radius
    tau = 9.6 * F * np.cbrt(F**2 * (2 - F)**2 * geo['d'] / Dm)
    side = np.pi * Dm / (geo['c'] * np.sin((alpha)*np.pi/180)) - (geo['n']) * geo['d'] / np.cos((90 - 2 * (alpha))*np.pi/180)
    l = 2 * side * np.cos((alpha)*np.pi/180)
    w = 2 * side * np.sin((alpha)*np.pi/180)
    e = np.sqrt(1 - (w/l)**2)
    vm = (2 / (3*np.sqrt(np.pi))) * (l/w)**(3/2) * ((1 - e**2) * e**2) / (ellipe(e) - (1 - e**2) * ellipk(e))
    Sr = l * w / 2
    m = Sr**(3/2) * vm
    Mh = 0.875 * u * m / (2 * np.pi**2 * Dm**2) * np.exp(-tau)
    
    nu = 2 * np.pi * Dm * np.sin((alpha)*np.pi/180) * np.cos((alpha)*np.pi/180) * F**2 / ((geo['n'])**2 * geo['d']**2)  #holes per unit length
    Lh = nu * Mh
    
    b = geo['n'] * geo['d'] * (1 - F) / F
    h = 2 * geo['d']**2 / (b + geo['d'])
    
    print(f'Shield coverage (K): {round(K * 100, 1)}%')
    print(f'Carrier separation (h): {h}mm ({round(h/geo["d"], 2)}d)')
    print(f'Hole inductance (Mh): {Mh}')
    print(f'Total inductance (Lh): {Lh}')
    print('Hole width: ', (W - geo['n']*geo['d']) / np.cos((alpha) * np.pi/180))
    
    print(W, F)
    
    return K, Mh, Lh

K, Mh, Lh = calc_params(geo)
