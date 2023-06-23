# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:55:30 2023

@author: griffin.kowash
"""

import numpy as np

class Vector:
    @staticmethod
    def mag(a):
        if a.ndim == 1:
            return np.sqrt(np.sum(a ** 2))
        else:
            return np.sqrt(np.sum(a ** 2, axis=1))
    
    @staticmethod
    def norm(a):
        if a.ndim == 1:
            return a / Vector.mag(a)
        else:
            return a / Vector.mag(a)[:, np.newaxis]