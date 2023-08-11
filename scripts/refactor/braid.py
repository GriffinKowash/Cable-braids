# -*- coding: utf-8 -*-
"""
Created on Fri May 26 08:49:21 2023

@author: griffin.kowash
"""

import os, time, json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from path import Path
from spline import Spline

            
class Braid:
    default_geo = {
        's': 1.00,
        'wave_width': 0,
        'wave_sharpness': 1,
        'c': 16,
        'n': 4,
        'd': 0.131,
        'alpha_deg': 72,
        'twist_deg': 7.53,
        'gap_shift_frac': 1/4,
        'tanh_c': 0.87,
        'tanh_m': 6.3
        }
    
    def __init__(self):
        self.geo = Braid.default_geo
        
        self.path = Path()
            
        self.curves = []
        self.data = []
                
    def set_geometry(self, geo=None, **kwargs):
        if geo != None:
            self.geo.update(geo)
            
        self.geo.update(kwargs)
                
        for key in self.geo.keys():
            setattr(self, key, self.geo[key])
            
        self.set_derived_params()
            
    def set_derived_params(self):
        self.alpha = self.alpha_deg * np.pi / 180
        self.twist = self.twist_deg * np.pi / 180
        self.w = 1 / (self.s * np.tan(self.alpha))                   # angular spatial frequency of helix (radians/mm)
        self.p = self.c / (4 * np.pi * np.tan(self.alpha) * self.s)  # "picks"/crossing frequency (1/mm)
        self.pitch = 2 * np.pi * self.s * np.tan(self.alpha)         # pitch height (mm)
        self.carrier_width = self.d * self.n                         # width of each carrier (mm)
        self.wave_w = 2 * np.pi * self.p / 2                         # angular frequency of wire waves with respect to length (radians / mm)
    
    def set_path_from_spline_file(self, filepath, resolution, equidistant=False): #deprecated
        spline = Spline()
        spline.set_from_file(filepath)
        points = spline.get_points(resolution)
        self.path.set_from_points(points, equidistant)
        
    def set_path_from_spline(self, params, resolution, equidistant=False):
        spline = Spline(params)
        points = spline.get_points(resolution)
        self.path.set_from_points(points, equidistant)
    
    def set_path_from_function(self, func, t_range, resolution, equidistant=False):
        self.path.set_from_function(func, t_range, resolution, equidistant)
    
    def set_path_from_points(self, points, equidistant=False):
        self.path.set_from_points(points, equidistant)
    
    def set_path_from_file(self, filepath, equidistant=False):
        self.path.set_from_file()
        
    def set_linear_path_between(self, p0, p1, resolution):
        self.path.set_linear_between(p0, p1, resolution)
    
    def construct(self, verbose=False, testing=False):     
        carrier_phis = np.linspace(self.twist, self.twist + 2*np.pi, int(self.c/2), endpoint=False)
        apply_power = lambda xs, p: np.array([x**p if x >= 0 else -(-x)**p for x in xs]) 
        
        def weave_func(ts, c=0.87, m=6.3):
            ts = ts % (2*np.pi)
            vals = []
            for t in ts:
                if t < c:
                    val = np.tanh(m * t) / np.tanh(m * c)
                elif c <= t < np.pi - c:
                    val = 1
                elif np.pi - c <= t < np.pi + c:
                    val = -np.tanh(m * (t - np.pi)) / np.tanh(m * c)
                elif np.pi + c <= t < 2*np.pi - c:
                    val = -1
                elif t > 2*np.pi - c:
                    val = np.tanh(m * (t - 2*np.pi)) / np.tanh(m * c)
                else:
                    print('Invalid value of t mod 2pi in weave_func: ', t)
                vals.append(val)
                
            return np.array(vals)
        
        gap_shift = 2*np.pi * self.gap_shift_frac * (np.pi * self.s / (self.c / 2 * np.cos(self.alpha)) - self.carrier_width) * np.sin(self.alpha) / (2*np.pi / self.wave_w)
        print(gap_shift)
        
        # loop over cw and ccw directions
        for chirality, sign in (('cw', -1), ('ccw', 1)):
            if verbose:
                print(f'\n\n_____ {chirality} carriers _____')
            
            # loop over carriers
            for i, carrier_phi in enumerate(carrier_phis):
                if verbose:
                    print(f'\n\tCARRIER {i}')
                
                phi_offset = (self.d / np.sin(self.alpha)) * (self.n - 1) / (2 * self.s)
                wire_phis = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, self.n, endpoint=True)
                
                if sign == 1:
                    if i % 2 == 0:
                        wave_offset = np.pi / 4
                    elif i % 2 == 1:
                        wave_offset = 5 * np.pi / 4
                elif sign == -1:
                    if i % 2 == 0:
                        wave_offset = -np.pi / 4
                    elif i % 2 == 1:
                        wave_offset = -5 * np.pi / 4
                
                # loop over wires
                for j, wire_phi in enumerate(wire_phis):
                    if verbose:
                        print(f'\t\twire {j}')
                    
                    # calculate position data
                    
                    #wire_wave_offset = (j - (self.n - 1) / 2) * (self.d / self.s) * np.sin(self.alpha) * np.cos(self.alpha)
                    wire_wave_offset = sign * (j - (self.n - 1) / 2) * (2 * np.pi * self.d * self.p)
                    
                    if testing:
                        wire_wave_offset += np.random.uniform(-np.pi, np.pi, 1)[0]
                    
                    #print(wire_wave_offset)
                    
                    #radius_offset = self.wave_width * apply_power(np.sin(self.wave_w * self.path.node_lengths + wave_offset + wire_wave_offset), self.wave_sharpness)

                    radius_offset = self.wave_width * weave_func(self.wave_w * self.path.node_lengths + wave_offset + wire_wave_offset + gap_shift, self.tanh_c, self.tanh_m)

                    s = self.s + radius_offset
                    x = s * np.cos(sign * self.w * self.path.node_lengths + wire_phi)
                    y = s * np.sin(sign * self.w * self.path.node_lengths + wire_phi)
                    r = np.array([x, y, np.zeros(self.path.resolution)])

                    # apply rotations and fit to path
                    for k, point in enumerate(r.T):
                        rotation = self.path.rotations[k]
                        new_point = self.path.nodes[k] + np.dot(rotation, np.array([x[k], y[k], 0]))
                        r.T[k] = new_point
                    
                    # store full data and nicely formatted curves
                    self.data.append((chirality, i, j, r.T))
                    self.curves.append((r[0,:], r[1,:], r[2,:]))
        
    def save(self, output_dir, fmt='default'):
        Saving.save(self, output_dir, fmt)
        
    def plot(self, mode='lines', linewidth=1.5):
        Plotting.plot(self, mode, linewidth)


class Plotting:
    @classmethod
    def plot(cls, braid, mode='lines', linewidth=1.5, ortho=True):
        fig = plt.figure()
        axes = fig.add_subplot(projection='3d')
        
        mins, maxs = cls.get_range(braid)
        
        if mode == 'lines':
            cls.plot_curves(braid, axes, linewidth)
            
        elif mode == 'surfaces':
            for carrier in range(braid.c):
                r0 = braid.data[carrier * braid.n][3]
                r1 = braid.data[carrier * braid.n + braid.n - 1][3]
                cls.plot_surface_between(axes, r0, r1)
        
        axes.set_xlim(mins[0], maxs[0])
        axes.set_ylim(mins[1], maxs[1])
        #axes.set_xlim(-1.2, 1.2)
        #axes.set_ylim(-1.2, 1.2)
        axes.set_zlim(mins[2], maxs[2])
        axes.view_init(elev=90, azim=0, roll=0)
        if ortho:
            axes.set_proj_type('ortho')
                
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
            
    @staticmethod
    def plot_surface_between(axes, r0, r1):
        for i in range(len(r0) - 1):
            vertices = [r0[i,:], r0[i+1,:], r1[i+1,:], r1[i,:]]
            axes.add_collection3d(Poly3DCollection([vertices], color='C0', alpha=0.6))
            

class Saving:
    @classmethod
    def save(cls, braid, output_dir, fmt='default'):
        output_dir = cls.check_output_dir(output_dir)
        
        if fmt == 'json':
            config = cls.get_config(braid, fmt)
            data = {'config': config, 'data': braid.data}
            with open(f'{output_dir}\\braid.json', 'w') as outfile:
                json.dump(data, outfile, cls=NumpyArrayEncoder)
        
        elif fmt == 'default':
            cls.save_config(braid, output_dir)
            for chirality, carrier, wire, data in braid.data:
                np.savetxt(f'{output_dir}\\data_{chirality}_c{carrier}_w{wire}.csv', data, delimiter=',')
    
    @staticmethod
    def check_output_dir(output_dir):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            print(f'Created output directory at {output_dir}')
            
        else:
            if len(os.listdir(output_dir)) > 2:  # this check is broken for save files in JSON format
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
                
        return output_dir
        
    @staticmethod
    def get_config(braid, fmt='default'):
        if fmt == 'default':
            return f"""carriers,{braid.c}
wires per carrier,{braid.n}
shield radius (mm),{braid.s}
wire diameter (mm),{braid.d}
carrier width (mm),{braid.carrier_width}
crossing frequency (1/mm),{braid.p}
pitch angle (degrees),{braid.alpha * 180 / np.pi}
pitch height (mm),{braid.pitch}
points per helix,{braid.path.resolution}"""
    
        elif fmt == 'json':
            return [('carriers', braid.c),
                    ('wires per carrier', braid.n),
                    ('shield radius (mm)', braid.s),
                    ('wire diameter (mm)', braid.d),
                    ('carrier width (mm)', braid.carrier_width),
                    ('crossing frequency (1/mm)', braid.p),
                    ('pitch angle (degrees)', braid.alpha_deg),
                    ('pitch height (mm)', braid.pitch),
                    ('points per helix', braid.path.resolution)]
        
    @classmethod
    def save_config(cls, braid, output_dir):
        config_text = cls.get_config(braid)
        config_file = open(output_dir + '\\config.txt', 'w')
        config_file.write(config_text)
        config_file.close()
          

class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
        