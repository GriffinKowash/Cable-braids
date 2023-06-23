# -*- coding: utf-8 -*-
"""
Created on Fri May 26 08:49:21 2023

@author: griffin.kowash
"""

import os, time
import numpy as np
import matplotlib.pyplot as plt
from path import Path

            
class Braid:
    default_geo = {
        's': 1.78/2,
        'c': 16,
        'n': 4,
        'd': 0.131,
        'alpha_deg': 72,
        }
    
    def __init__(self):
        self.geo = Braid.default_geo
        
        self.path = Path()
            
        self.curves = []
        self.surfaces = []
        self.data = []
                
    def set_geometry(self, geo=None, **kwargs):
        if geo != None:
            self.geo.update(geo)
            
        self.geo.update(kwargs)
                
        for key in self.geo.keys():
            setattr(self, key, self.geo[key])
            
        self.set_derived_params()
            
    def set_derived_params(self):
        self.alpha = self.geo['alpha_deg'] * np.pi / 180
        self.w = 1 / (self.s * np.tan(self.alpha))                   # angular frequency with respect to length (radians/mm)
        self.p = np.tan(self.alpha) * self.c / (4 * np.pi * self.s)  # "picks"/crossing frequency (1/mm)
        self.pitch = 8 * np.pi**2 *  self.s**2 * self.p / self.c     # pitch height (mm)
        self.carrier_width = self.d * self.n                         # width of each carrier (mm)
    
    def set_path_from_spline(self):
        pass
    
    def set_path_from_function(self, func, t_range, resolution, equidistant=False):
        self.path.set_from_function(func, t_range, resolution, equidistant)
    
    def set_path_from_points(self, points, equidistant=False):
        self.path.set_from_points(points, equidistant)
    
    def set_path_from_file(self, filepath, equidistant=False):
        self.path.set_from_file()
        
    def set_linear_path_between(self, p0, p1, resolution):
        self.path.set_linear_between(p0, p1, resolution)
    
    
    def construct(self, verbose=False):
        # Calculate braid
        t = np.linspace(0, self.path.length, self.path.resolution, endpoint=True)
        z = t
        carrier_phis = np.linspace(0, 2*np.pi, int(self.c/2), endpoint=False)
        
        # loop over cw and ccw directions
        for chirality, sign in (('cw', -1), ('ccw', 1)):
            if verbose:
                print(f'\n\n_____ {chirality} carriers _____')
            
            # loop over carriers
            for i, carrier_phi in enumerate(carrier_phis):
                if verbose:
                    print(f'\n\tCARRIER {i}')
                
                phi_offset = self.d * (self.n - 1) / (2 * self.s)  #should include np.sin(self.alpha) term in denominator?
                wire_phis = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, self.n, endpoint=True)
                
                # loop over wires
                for j, wire_phi in enumerate(wire_phis):
                    if verbose:
                        print(f'\t\twire {j}')
                    
                    # calculate position data
                    x = self.s * np.cos(sign * self.w * t + wire_phi)
                    y = self.s * np.sin(sign * self.w * t + wire_phi)
                    r = np.array([x, y, z])
                                        
                    # apply rotation if following a path
                    if self.path != None:
                        for k, point in enumerate(r.T):
                            rotation = self.path.rotations[k]
                            new_point = self.path.nodes[k] + np.dot(rotation, np.array([x[k], y[k], 0]))
                            r.T[k] = new_point
                    
                    # store save data
                    self.data.append((chirality, i, j, r.T))
    
                    # store curve data (could derive from self.data)
                    self.curves.append((r[0,:], r[1,:], r[2,:]))
                
                # store surface data (could derive from self.data)
                x_grid, y_grid, z_grid = self.generate_surface(carrier_phi, phi_offset, sign)
                self.surfaces.append((x_grid, y_grid, z_grid))
        
    def generate_surface(self, carrier_phi, phi_offset, sign):
        z_vals = np.linspace(0, self.path.length, self.path.resolution)
        phi_vals = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, 10)
        phi_grid, z_grid = np.meshgrid(phi_vals, z_vals)
        
        phi_shift = np.reshape(self.w * z_vals, (-1,1))  # shift in phi as z-position increases
        phi_grid = phi_grid + phi_shift
        
        x_grid = self.s * np.cos(sign * phi_grid)
        y_grid = self.s * np.sin(sign * phi_grid)
        
        return x_grid, y_grid, z_grid
    
    def save(self, output_dir):
        Saving.save(self, output_dir)
    
    def plot(self, mode='lines', linewidth=1.5):
        Plotting.plot(self, mode, linewidth)


class Plotting:
    @classmethod
    def plot(cls, braid, mode='lines', linewidth=1.5):
        fig = plt.figure()
        axes = fig.add_subplot(projection='3d')
        
        mins, maxs = cls.get_range(braid)
        
        if mode == 'lines':
            cls.plot_curves(braid, axes, linewidth)
            
        elif mode == 'surfaces':
            cls.plot_surfaces(braid, axes)
        
        axes.set_xlim(mins[0], maxs[0])
        axes.set_ylim(mins[1], maxs[1])
        axes.set_zlim(mins[2], maxs[2])
        
        fig.show()
        
    @staticmethod
    def get_range(braid, buffer=0.05, square=True):
        curves = np.array(braid.curves)
        points = curves.swapaxes(1, 2).reshape(curves.shape[0] * curves.shape[2], 3)
        
        mid = points.mean(axis=0)
        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        deltas = (maxs - mins) * (1 + buffer) / 2
        
        if square:
            deltas = np.max(deltas) * np.ones(3)
        
        mins = mid - deltas
        maxs = mid + deltas
                
        return mins, maxs
        
    @staticmethod
    def plot_curves(braid, axes, linewidth):
        for x, y, z in braid.curves:
            axes.plot(x, y, z, linewidth=linewidth, color='C0')
            
    @staticmethod
    def plot_surfaces(braid, axes):
        for x_grid, y_grid, z_grid in braid.surfaces:
            axes.plot_surface(x_grid, y_grid, z_grid, color=(0.12, 0.47, 0.71, 0.4))
            

class Saving:
    @classmethod
    def save(cls, braid, output_dir):
        cls.check_output_dir(output_dir)
        cls.save_config(braid, output_dir)
        
        for chirality, carrier, wire, data in braid.data:
            np.savetxt(f'{braid.params.output_dir}\\data_{chirality}_c{carrier}_w{wire}.csv', data, delimiter=',')
    
    @staticmethod
    def check_output_dir(output_dir):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            print(f'Created output directory at {output_dir}')
            
        else:
            if len(os.listdir(output_dir)) > 2:
                if 'y' in input('Output directory is not empty. Overwrite data? [y/n] ').lower():
                    for file_name in os.listdir(output_dir):
                        file = output_dir + file_name
                        if os.path.isfile(file):
                            os.remove(file)
                    print(f'Overwriting data at {output_dir}')
                    
                else:
                    output_dir = output_dir + '_' + str(int(time.time()))
                    os.mkdir(output_dir)
                    print(f'Created new output directory at {output_dir}')
                    
            else:
                print(f'Outputting data to directory {output_dir}')
        
    @staticmethod
    def save_config(braid, output_dir):
        config_text = \
f"""carriers,{braid.c}
wires per carrier,{braid.n}
shield radius (mm),{braid.s}
wire diameter (mm),{braid.d}
carrier width (mm),{braid.carrier_width}
crossing frequency (1/mm),{braid.p}
pitch angle (degrees),{braid.alpha * 180 / np.pi}
pitch height (mm),{braid.pitch}
points per helix,{braid.resolution}"""

        config_file = open(output_dir + '\\config.txt', 'w')
        config_file.write(config_text)
        config_file.close()
          

        
        