# -*- coding: utf-8 -*-
"""
Created on Fri May 26 08:49:21 2023

@author: griffin.kowash
"""

import os, time, json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from path import Path
from spline import Spline

            
class Braid:
    default_geo = {
        'mode': 'body',      # modeling style (line, surface, body)
        'ID': 1.00,             # inner diameter (mm)
        'c': 16,                # number of carriers
        'n': 4,                 # wires per carrier
        'd': 0.127,             # wire diameter (mm)
        'h': None,              # carrier separation ("None" for auto-calculate)
        'alpha_deg': 18,        # weave angle (degrees; relative to cable axis)
        'twist_deg': 0,         # overall rotation (degrees; seldom modified from default)
        'weave_offset': 0,      # allows offsetting of weave function
        'tanh_c': 0.5,          # weave function parameter
        'tanh_m': 1,            # weave function parameter
        'wire_spread': 1,       # artificial wire separation (fraction typically 1.0~1.1; only for body models)
        'skip_factor': 1        # affects carrier resolution in surface models
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
        if self.h == None:
            self.h = self.calc_h()                                              # height between carriers as defined by Tyni
        self.Dm = self.calc_Dm()                                                # mean diameter (center of carrier separation)
        self.w = 1 / (self.Dm / 2 * (1 / np.tan(self.alpha)))                   # angular spatial frequency of helix (radians/mm)
        self.p = self.c / (2 * np.pi * (1 / np.tan(self.alpha)) * self.Dm)      # "picks"/crossing frequency (1/mm)
        self.pitch = np.pi * self.Dm * (1 / np.tan(self.alpha))                 # pitch height (mm)
        self.carrier_width = self.d * self.n                                    # width of each carrier (mm)
        self.weave_w = 2 * np.pi * self.p / 2                                   # angular frequency of weave pattern with respect to length (radians / mm)
        
        if self.mode == 'body':                                                 # ampltiude of weave pattern (mean radius as reference)
            self.amplitude = (self.d + self.h) / 2                              # (extend to centerlines for body models)
        elif self.mode == 'surface' or self.mode == 'line':                     # (extend to inner edge for line/surface models)
            self.amplitude = (self.d + self.h) / 2                              
        else:
            print('Unrecognized mode ', self.mode)

    def calc_Dm(self):
        return self.ID + 2*self.d + self.h

    def calc_h(self):
        W = 2 * np.pi * self.ID / self.c * np.cos(self.alpha)
        F = self.n * self.d / W
        b = self.n * self.d * (1 - F) / F
        h = 2 * self.d**2 / (b + self.d)
        print('h = ', h)
        return h
    
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
        
        gap_shift = 2*np.pi * self.weave_offset * (np.pi * self.Dm / (self.c * np.sin(self.alpha)) - self.carrier_width) * np.cos(self.alpha) / (2*np.pi / self.weave_w)  #aligns up/down weaving with gaps between carriers (ideally)
        
        
        # loop over cw and ccw directions
        for chirality, sign in (('cw', -1), ('ccw', 1)):
            if verbose:
                print(f'\n\n_____ {chirality} carriers _____')
            
            # loop over carriers
            for i, carrier_phi in enumerate(carrier_phis):
                if verbose:
                    print(f'\n\tCARRIER {i}')
                    
                # calculate angular positions of wires within carrier
                if self.mode == 'body':
                    # body: wires modeled along centerlines
                    n_lines = self.n
                    phi_offset = (self.d / np.cos(self.alpha)) * self.n / self.Dm * self.wire_spread
                    wire_phis = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, self.n, endpoint=True)
                
                elif self.mode == 'line':
                    # line: wires modeled out to carrier edges to ensure correct aperture size
                    n_lines = self.n
                    phi_offset = (self.d / np.cos(self.alpha)) * (self.n + 1) / self.Dm
                    
                elif self.mode == 'surface':
                    # surface: wires modeled out to carrier edges with optional skip factor
                    n_lines = round((self.n + 1) / self.skip_factor)
                    phi_offset = (self.d / np.cos(self.alpha)) * (self.n + 1) / self.Dm
                    
                wire_phis = np.linspace(carrier_phi - phi_offset, carrier_phi + phi_offset, n_lines, endpoint=True)
        
                # calculate phase offset for each carrier to avoid intersection
                if sign == 1:
                    if i % 2 == 0:
                        carrier_weave_offset = np.pi / 4
                    elif i % 2 == 1:
                        carrier_weave_offset = 5 * np.pi / 4
                elif sign == -1:
                    if i % 2 == 0:
                        carrier_weave_offset = -np.pi / 4
                    elif i % 2 == 1:
                        carrier_weave_offset = -5 * np.pi / 4
                
                # loop over wires
                for j, wire_phi in enumerate(wire_phis):
                    if verbose:
                        print(f'\t\twire {j}')
                    
                    # calculate weave phase offset for each wire to align inflection points with gaps between carriers
                    #wire_weave_offset = -sign * (j - (n_lines - 1) / 2) * (self.weave_w * self.d * np.sin(self.alpha) / np.tan(np.pi - 2*self.alpha))  #double check change of alpha definition
                    wire_weave_offset = sign * self.Dm / (4 * np.tan(self.alpha)) * (wire_phi - carrier_phi) * self.weave_w
                    
                    #if self.mode == 'line':
                    #    wire_weave_offset = -sign * (j - (self.n - 1) / 2) * (self.weave_w * self.d * np.sin(self.alpha) / np.tan(np.pi - 2*self.alpha))  #double check change of alpha definition
                    #    
                    #elif self.mode == 'surface':
                    #    if j == 0:
                    #        wire_weave_offset = sign * (0 - (self.n - 1) / 2) * (self.d * self.weave_w / np.cos(2 * self.alpha))
                    #    elif j == 1:
                    #        wire_weave_offset = sign * ((self.n - 1) - (self.n - 1) / 2) * (self.d * self.weave_w / np.cos(2 * self.alpha))
                       
                    radius_offset = self.amplitude * weave_func(self.weave_w * self.path.node_lengths + carrier_weave_offset + wire_weave_offset + gap_shift, self.tanh_c, self.tanh_m)

                    s = self.Dm / 2 + radius_offset
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
        
        
    def save(self, output_dir, fmt='default', mode='line'):
        Saving.save(self, output_dir, fmt, mode)
        
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
    def save(cls, braid, output_dir, fmt='default', mode='line'):
        output_dir = cls.check_output_dir(output_dir)
        
        if fmt == 'json':
            config = cls.get_config(braid, fmt, mode)
            data = {'config': config, 'data': braid.data}
            with open(f'{output_dir}\\braid.json', 'w') as outfile:
                json.dump(data, outfile, cls=NumpyArrayEncoder)
        
        elif fmt == 'default':
            cls.save_config(braid, output_dir, mode)
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
    def get_config(braid, fmt='default', mode='line'):
        if fmt == 'default':
            return f"""carriers,{braid.c}
wires per carrier,{braid.n}
inner diameter (mm),{braid.ID}
mean diameter (mm),{braid.Dm}
wire diameter (mm),{braid.d}
weave angle (degrees),{braid.alpha * 180 / np.pi}
carrier separation (mm):,{braid.h}
points per helix,{braid.path.resolution}"""
    
        elif fmt == 'json':
            return [('carriers', braid.c),
                    ('wires per carrier', braid.n),
                    ('inner diameter (mm)', braid.ID),
                    ('mean diameter (mm)', braid.Dm),
                    ('wire diameter (mm)', braid.d),
                    ('weave angle (degrees)', braid.alpha_deg),
                    ('carrier separation (mm)', braid.h),
                    ('points per helix', braid.path.resolution)]
        
    @classmethod
    def save_config(cls, braid, output_dir, mode='line'):
        config_text = cls.get_config(braid, fmt='default', mode=mode)
        config_file = open(output_dir + '\\config.txt', 'w')
        config_file.write(config_text)
        config_file.close()
          

class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
        