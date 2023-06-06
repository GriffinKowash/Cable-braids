# -*- coding: utf-8 -*-
"""
Created on Fri May 26 08:49:21 2023

@author: griffin.kowash
"""

import sys, os, time
import getopt
import numpy as np
import matplotlib.pyplot as plt

class Params:
    # Constants
    PICKS = 0
    PITCH_ANGLE = 1
    LINES = 0
    SURFACES = 1
    
    def __init__(self):
        # Settings
        self.output_dir = r'C:\Users\griffin.kowash\AppData\Local\Temp\braids\braid_data'
        self.param_mode = Params.PITCH_ANGLE  # determine parameters using either pitch angle or picks (crossings per mm)
        self.plot_mode = Params.LINES         # plot showing individual wires (LINES) or carrier ribbons (SURFACES)
        self.plotting = False                 # visualize braid in matplotlib
        self.saving = True                   # save output files for use in Discovery
        self.verbose = False                  # provide detailed output to monitor braid construction
    
        # Input parameters (only alpha *or* p required; choose using "param_mode" above)
        self.s = 1.78 / 2       # shield radius (mm)
        self.c = 16             # number of carriers (count)
        self.n = 4              # "ends"/wires per carrier (count)
        self.d = 0.131          # wire diameter (mm)
        self.p = 100            # "picks"/crossing frequency (1/mm)
        self.alpha = 72         # pitch angle (degrees)
        self.z_max = 145        # endpoint of z axis (mm)
        self.resolution = 150   # number of points per curve (count)

        # Parse command line arguments
        if len(sys.argv) > 1:
            self.set_from_args()
            
        # Set path along which to construct braid (default along z axis)
        path_func = lambda t: (0.00001*t, 0.00001*t, t)
        self.set_path(path_func, (0, self.z_max))
        
        # Calculate derived paramaters
        self.get_derived_params()
        
    def set_from_args(self):
        # Command line parameters
        arg_str = "s:c:n:d:a:z:r:p:m:f:v:o:"
        arg_list = ['plot', 'plotmode =', 'save', 'verbose', 'outdir =']
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
            else:
                print(f'Unknown option {opt} provided.')
    
    def get_derived_params(self):
        if self.param_mode == Params.PICKS:
            self.alpha = np.arctan(4*np.pi * self.s * self.p / self.c)   # pitch angle (radians)
        elif self.param_mode == Params.PITCH_ANGLE:
            self.alpha *= np.pi / 180                                    # pitch angle (radians)
            self.p = np.tan(self.alpha) * self.c / (4 * np.pi * self.s)  # "picks"/crossing frequency (1/mm)
            
        self.pitch = 8 * np.pi**2 *  self.s**2 * self.p / self.c         # pitch height (mm)
        self.w = self.c / (4 * np.pi * self.s**2 * self.p)               # angular frequency with respect to z (radians/mm)
        self.carrier_width = self.d * self.n                             # width of each carrier (mm)

    def restrict_params(self):
        if self.c % 2 == 1:
            print('Odd number of carriers provided for parameter c. Assuming that value refers to number of carriers per direction.')
            self.c *= 2
            
    def set_path(self, func, t_range):
        self.path = Path()
        self.path.set_from_function(func, t_range, self.resolution + 1)  ### this could cause issues if Params.resolution is changed after creating the path.
            
            
class Braid:
    def __init__(self, params):
        self.params = params
        
        self.curves = []
        self.surfaces = []
        self.data = []
        
        self.construct()
        
        if self.params.plotting:
            self.plot()
            
        if self.params.saving:
            self.save()
            
    def construct(self):
        # Calculate braid, create plot, and save data
        t = np.linspace(0, self.params.z_max, self.params.resolution)
        z = t
        carrier_phis = np.linspace(0, 2*np.pi, int(self.params.c/2), endpoint=False)
        
        # loop over cw and ccw directions
        for chirality, sign in (('cw', -1), ('ccw', 1)):
            if self.params.verbose:
                print(f'\n\n_____ {chirality} carriers _____')
            
            # loop over carriers
            for i, carrier_phi in enumerate(carrier_phis):
                if self.params.verbose:
                    print(f'\n\tCARRIER {i}')
                
                phi_offset = self.params.d * (self.params.n - 1) / (2 * self.params.s)  #should include np.sin(self.params.alpha) term in denominator
                wire_phis = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, self.params.n, endpoint=True)
                
                # loop over wires
                for j, wire_phi in enumerate(wire_phis):
                    if self.params.verbose:
                        print(f'\t\twire {j}')
                    
                    # calculate position data
                    x = self.params.s * np.cos(sign * self.params.w * t + wire_phi)
                    y = self.params.s * np.sin(sign * self.params.w * t + wire_phi)
                    r = np.array([x, y, z])
                    
                    # apply rotation if following curve
                    if self.params.path != None:
                        for k, point in enumerate(r.T):
                            rotation = self.params.path.rotations[k]
                            new_point = self.params.path.nodes[k] + np.dot(rotation, np.array([x[k], y[k], 0]))
                            r.T[k] = new_point
                    
                    
                    # store save data
                    self.data.append((chirality, i, j, r.T))
    
                    # store curve data (could derive from self.data)
                    self.curves.append((r[0,:], r[1,:], r[2,:]))
                
                # store surface data (could derive from self.data)
                x_grid, y_grid, z_grid = self.generate_surface(carrier_phi, phi_offset, sign)
                self.surfaces.append((x_grid, y_grid, z_grid))
        
    def generate_surface(self, carrier_phi, phi_offset, sign):
        z_vals = np.linspace(0, self.params.z_max, self.params.resolution)
        phi_vals = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, 10)
        phi_grid, z_grid = np.meshgrid(phi_vals, z_vals)
        
        phi_shift = np.reshape(self.params.w * z_vals, (-1,1))  # shift in phi as z-position increases
        phi_grid = phi_grid + phi_shift
        
        x_grid = self.params.s * np.cos(sign * phi_grid)
        y_grid = self.params.s * np.sin(sign * phi_grid)
        
        return x_grid, y_grid, z_grid
    
    def save(self):
        print('Calling save method of Braid.')
        Saving.save(self)
                
    def plot(self):
        Plotting.plot(self)
        
        
class Path:
    def __init__(self):
        pass
    
    def set_from_points(self, points):
       
        self.xhat = np.array([1,0,0])
        self.yhat = np.array([0,1,0])
        self.zhat = np.array([0,0,1])

        self.nodes = []
        self.normals = []
        self.rotations = []
        
        for p0, p1 in zip(points[:-1], points[1:]):
            node = (p0 + p1) / 2
            normal = Vector.norm(p1 - p0)
            rotation = self.get_rotation_matrix(self.zhat, normal)
            self.nodes.append(node)
            self.normals.append(normal)
            self.rotations.append(rotation)
    
    def set_from_function(self, func, t_range, resolution):
        t = np.linspace(t_range[0], t_range[1], resolution)
        points = np.array(func(t)).T
        self.set_from_points(points)
    
    def get_rotation_matrix(self, a, b):
        # Rotates z direction to align with helix axis
        G = np.zeros((3,3))
        G[0,0] = np.dot(a,b)
        G[1,0] = Vector.mag(np.cross(a, b))
        G[0,1] = -1*G[1,0]
        G[1,1] = G[0,0]
        G[2,2] = 1
        
        Fi = np.array([a, Vector.norm(b - np.dot(a, b) * a), np.cross(b, a)]).T
        F = np.linalg.inv(Fi)
        U = np.dot(Fi, np.dot(G, F))
        return U


class Plotting:
    @classmethod
    def plot(cls, braid):
        fig = plt.figure()
        axes = fig.add_subplot(projection='3d')
        
        if braid.params.plot_mode == Params.LINES:
            if braid.params.verbose:
                print('\nPlotting curves.')
            cls.plot_curves(braid, axes)
            
        elif braid.params.plot_mode == Params.SURFACES:
            if braid.params.verbose:
                print('\nPlotting surfaces.')
            cls.plot_surfaces(braid, axes)
            
        axes.set_xlim(-braid.params.z_max/2, braid.params.z_max/2)
        axes.set_ylim(-braid.params.z_max/2, braid.params.z_max/2)
        axes.set_zlim(0, braid.params.z_max)
        fig.show()
            
    @staticmethod
    def plot_curves(braid, axes):
        for x, y, z in braid.curves:
            axes.plot(x, y, z, color='C0')
            
    @staticmethod
    def plot_surfaces(braid, axes):
        for x_grid, y_grid, z_grid in braid.surfaces:
            axes.plot_surface(x_grid, y_grid, z_grid, color=(0.12, 0.47, 0.71, 0.4))
            

class Saving:
    @classmethod
    def save(cls, braid):
        cls.check_output_dir(braid)
        cls.save_config(braid)
        
        for chirality, carrier, wire, data in braid.data:
            np.savetxt(f'{braid.params.output_dir}\\data_{chirality}_c{carrier}_w{wire}.csv', data, delimiter=',')
    
    @staticmethod
    def check_output_dir(braid):
        if not os.path.exists(braid.params.output_dir):
            os.mkdir(braid.params.output_dir)
            print(f'Created output directory at {braid.params.output_dir}')
            
        else:
            if len(os.listdir(braid.params.output_dir)) != 0:
                if 'y' in input('Output directory is not empty. Overwrite data? [y/n] ').lower():
                    for file_name in os.listdir(braid.params.output_dir):
                        file = braid.params.output_dir + file_name
                        if os.path.isfile(file):
                            os.remove(file)
                    print(f'Overwriting data at {braid.params.output_dir}')
                    
                else:
                    braid.params.output_dir = braid.params.output_dir + '_' + str(int(time.time()))
                    os.mkdir(braid.params.output_dir)
                    print(f'Created new output directory at {braid.params.output_dir}')
                    
            else:
                print(f'Outputting data to directory {braid.params.output_dir}')
        
    @staticmethod
    def save_config(braid):
        config_text = \
f"""carriers,{braid.params.c}
wires per carrier,{braid.params.n}
shield radius (mm),{braid.params.s}
wire diameter (mm),{braid.params.d}
carrier width (mm),{braid.params.carrier_width}
crossing frequency (1/mm),{braid.params.p}
pitch angle (degrees),{braid.params.alpha * 180 / np.pi}
pitch height (mm),{braid.params.pitch}
points per helix,{braid.params.resolution}"""

        config_file = open(braid.params.output_dir + '\\config.txt', 'w')
        config_file.write(config_text)
        config_file.close()
          

class Vector:
    @staticmethod
    def mag(a):
        return np.sqrt(np.sum(a ** 2))
    
    @staticmethod
    def norm(a):
        return a / Vector.mag(a)

        
        
if __name__ == '__main__':
    path_func = lambda t: (4*16*np.sin(t)**3, 4*(13*np.cos(t) - 5*np.cos(2*t) - 2*np.cos(3*t) - np.cos(4*t)), 0*np.sin(t))
    #path = Path()
    #path.set_from_function(path_func, (0, 10), 151)
    
    params = Params()
    params.resolution = 2000
    params.set_path(path_func, (0, 2*np.pi))
    
    #params.verbose = True
    params.plotting = True
    params.saving = False
    #params.plot_mode = Params.SURFACES

    
    braid = Braid(params)