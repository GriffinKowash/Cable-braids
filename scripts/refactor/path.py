# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:54:36 2023

@author: griffin.kowash
"""

import numpy as np
import matplotlib.pyplot as plt
from vector import Vector

class Path:
    def __init__(self):
        self.func = None
    
    def set_linear_between(self, p0, p1, resolution):  
        p0, p1 = np.array(p0), np.array(p1)
        to_p1 = Vector.norm(p1 - p0)
        dist = Vector.mag(p1 - p0)
        
        points = []
        for t in np.linspace(0, 1, resolution):
            points.append(p0 + to_p1 * dist * t)
        
        self.set_from_points(points)
    
    def set_from_points(self, points, equidistant=False, testing=True):        
        if not testing:
            points = np.array(points)
            self.resolution = points.shape[0] - 1
            print('resolution: ', self.resolution)
            self.zhat = np.array([0,0,1])
    
            self.nodes = []
            self.normals = []
            self.rotations = []
                        
            for p0, p1 in zip(points[:-1], points[1:]):
                node = (p0 + p1) / 2
                #print(p0, p1)
                normal = Vector.norm(p1 - p0)
                #print(normal)
                rotation = self.get_rotation_matrix(self.zhat, normal)
                #print(rotation)
                self.nodes.append(node)
                self.normals.append(normal)
                self.rotations.append(rotation)
                
            self.nodes = np.array(self.nodes)
            self.normals = np.array(self.normals)
            self.rotations = np.array(self.rotations)
                
            self.length = self.get_total_length()
            self.node_lengths = self.get_node_lengths()
            
            if equidistant:
                self.calc_equidistant_nodes(points.shape[0])
                
        else:
            points = np.array(points)
            self.resolution = points.shape[0]
            print('resolution: ', self.resolution)
            self.zhat = np.array([0,0,1])
    
            self.nodes = []
            self.normals = []
            self.rotations = []
                        
            for i in range(self.resolution):                
                if i == 0:
                    p0, p1 = points[i], points[i+1]
                elif i == self.resolution - 1:
                    p0, p1 = points[i-1], points[i]     
                else:
                    p0, p1 = points[i-1], points[i+1]
                
                node = points[i]
                normal = Vector.norm(p1 - p0)  #note that "normal" is normal to desired helix plane; actually tangent to path. Confusing.
                rotation = self.get_rotation_matrix(self.zhat, normal)

                self.nodes.append(node)
                self.normals.append(normal)
                self.rotations.append(rotation)
                
            print('rotations: ', len(self.rotations))
                
            self.nodes = np.array(self.nodes)
            self.normals = np.array(self.normals)
            self.rotations = np.array(self.rotations)
                
            self.length = self.get_total_length()
            self.node_lengths = self.get_node_lengths()
            
            if equidistant:
                self.calc_equidistant_nodes()
    
    def set_from_function(self, func, t_range, resolution, equidistant=False):
        self.func = func
        self.t_range = t_range
                
        t = np.linspace(t_range[0], t_range[1], resolution)
        points = np.array(func(t)).T
        self.set_from_points(points, equidistant)
        
        
    def set_from_file(self, filepath, equidistant=False):
        points = np.loadtxt(filepath)
        self.set_from_points(points, equidistant)
    
    def get_rotation_matrix(self, a, b):
        # Rotates vector a to vector b
        if any(np.cross(a, b) != 0):
            G = np.zeros((3,3))
            G[0,0] = np.dot(a,b)
            G[1,0] = Vector.mag(np.cross(a, b))
            G[0,1] = -1*G[1,0]
            G[1,1] = G[0,0]
            G[2,2] = 1
            
            Fi = np.array([a, Vector.norm(b - np.dot(a, b) * a), np.cross(b, a)]).T
            F = np.linalg.inv(Fi)
            U = np.dot(Fi, np.dot(G, F))
        else:
            U = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
        return U
    
    def get_total_length(self):
        return self.get_length_between(0, -1)
    
    def get_length_between(self, t0, t1):
        if t1 == -1:
            diff = self.nodes[t0+1:, :] - self.nodes[t0:t1, :]
        else:
            #print(self.nodes.shape)
            #print(t0)
            #print(t1)
            diff = self.nodes[t0+1:t1+1, :] - self.nodes[t0:t1, :]
            
        return np.sum(np.sqrt(np.sum(np.power(diff, 2), axis=1)))
    
    def get_node_lengths(self):
        # Find distance from start point for each node
        node_lengths = []
        for i in range(len(self.nodes)):
            node_lengths.append(self.get_length_between(0, i))
        return node_lengths
            
    def refine(self, factor):
        # remake mesh with resolution increased by given integer factor
        if type(factor) != int:
            print(f'Non-integer refinement factor {factor} provided to Path.refine; setting value to {int(factor)}.')
            factor = int(factor)
            
        self.resolution *= factor
        
        print('nodes before refinement: ', self.nodes.shape)


        if self.func != None:
            self.set_from_function(self.func, self.t_range, self.resolution)  
                        
        else:
            new_nodes = [self.nodes[0]]
            for i in range(1, len(self.nodes)):
                n1, n0 = self.nodes[i], self.nodes[i-1]
                to_n1 = n1 - n0
                for j in range(1, factor + 1):
                    new_nodes.append(n0 + to_n1 * j / factor)
            self.nodes = np.array(new_nodes)
            
        self.length = self.get_total_length()
        self.node_lengths = self.get_node_lengths()
              
    def calc_equidistant_nodes(self):
        original_resolution = self.resolution
        
        interval = self.length / self.resolution
        
        largest_segment = max([Vector.mag(self.nodes[i] - self.nodes[i-1]) for i in range(1, len(self.nodes))])
        refine_factor = int(np.ceil(largest_segment / interval))
        
        self.nodes_old = self.nodes
        self.normals_old = self.normals
        self.rotations_old = self.rotations
        
        self.refine(refine_factor)
        
        print('refine ', refine_factor)
        
        self.resolution = original_resolution # refine method updates resolution, but since we only want refined path for calculating equidistant path, we reset it manually here
        interval = self.length / self.resolution
 
        self.nodes_eq = [self.nodes[0]]
        self.normals_eq = [self.normals[0]]
        self.rotations_eq = [self.rotations[0]]
        
        print('interval: ', interval)
        print('length: ', self.length)
        print('resolution: ', self.resolution)
        
        
        for i in range(1, self.resolution * refine_factor - 1): 
            target_length = len(self.nodes_eq) * interval
            length_since_start = self.get_length_between(0, i)
            
            if length_since_start >= target_length:
                n0, n1 = self.nodes[i-1], self.nodes[i]
                
                diff = length_since_start - target_length
                segment_length = self.get_length_between(i-1, i)
                frac = 1 - diff / segment_length
                
                new_node = (1 - frac) * n0 + frac * n1
                new_normal = Vector.norm(new_node - n0)
                new_rotation = self.get_rotation_matrix(self.zhat, new_normal)

                self.nodes_eq.append(new_node)
                self.normals_eq.append(new_normal)
                self.rotations_eq.append(new_rotation)

        self.nodes_ref = np.array(self.nodes)
        self.normals_ref = np.array(self.normals)
        self.rotations_ref = np.array(self.rotations)
                
        self.nodes = np.array(self.nodes_eq)
        self.normals = np.array(self.normals_eq)
        self.rotations = np.array(self.rotations_eq)
        
        self.length = self.get_total_length()
        self.node_lengths = self.get_node_lengths()
                
    def plot_xy_nodes(self, skip=1, equidistant=False):
        fig, ax = plt.subplots(1)
        
        if equidistant:
            ax.scatter(self.nodes_ref[::skip, 0], self.nodes_ref[::skip, 1], marker='.')
        ax.scatter(self.nodes[::skip, 0], self.nodes[::skip, 1], marker='.')

        fig.show()
        
    def plot_node_spacing(self, equidistant=False):
        fig, ax = plt.subplots(1)
        
        if equidistant:
           old_spacing = Vector.mag(self.nodes_old[1:] - self.nodes_old[:-1])
           ax.plot(range(len(self.nodes_old) - 1), old_spacing)
           
        spacing = Vector.mag(self.nodes[1:] - self.nodes[:-1])
        ax.plot(range(len(self.nodes) - 1), spacing)
        #ax.set_yscale('log')
        fig.show()
        
    def plot_node_spacing_deviation(self, equidistant=False):
        fig, ax = plt.subplots(1)
        
        if equidistant:
           old_spacing = Vector.mag(self.nodes_old[1:] - self.nodes_old[:-1])
           old_spacing_dev = old_spacing / old_spacing.mean() - 1
           ax.plot(range(len(self.nodes_old) - 1), old_spacing_dev)
           #print('Maximum deviation as fraction of mean (original): ', max(abs(old_spacing_dev)) * 100, '%')
           print('Standard deviation (original): ', np.std(old_spacing))
           
        spacing = Vector.mag(self.nodes[1:] - self.nodes[:-1])
        spacing_dev = spacing / spacing.mean() - 1
        #print('Maximum deviation as fraction of mean: ', max(abs(spacing_dev)) * 100, '%')
        print('Standard deviation (current): ', np.std(spacing))
        ax.plot(range(len(self.nodes) - 1), spacing_dev)
        #ax.set_yscale('log')
        fig.show()