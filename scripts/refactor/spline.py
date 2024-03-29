from geomdl import BSpline, convert, fitting, exchange
import numpy as np
import matplotlib.pyplot as plt


class Spline:
    def __init__(self, params=None):
        if params != None:
            for key in params.keys():
                setattr(self, key, params[key])
    
    def set_from_file(self, filepath):
        # to-do: update for json
        with open(filepath, 'r') as file:
            lines = file.readlines()

        self.degree = int(lines[0].split(': ')[1][:-1])
        self.knot_vector = [float(knot) for knot in lines[1].split(': ')[1][:-1].split(', ')]

        self.control_points = []
        for line in lines[2:]:
            x, y, z = line.split('\n')[0].split(', ')
            self.control_points.append((1000*float(x), 1000*float(y), 1000*float(z)))
        
    def set_from_dict(self, params):
        for key in params.keys():
            setattr(self, key, params[key])

    def get_points(self, resolution):
        crv = BSpline.Curve()
        crv.degree = self.degree
        crv.ctrlpts = self.control_points
        crv.knotvector = self.knot_vector

        crv_rat = convert.bspline_to_nurbs(crv)

        crv_rat.sample_size = resolution
        points = np.array(crv_rat.evalpts)
        
        return points
    
    def interpolate_points(self, points):
        #points = exchange.import_csv(points)
        curve = fitting.interpolate_curve(points, degree=3)
        self.control_points = curve.ctrlpts
        self.knot_vector = curve.knotvector
        self.degree = curve.degree
    
    def plot(self, fig=None, ax=None, show=True):
        if fig == None or ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
        
        points = self.get_points(300)
        control_points = np.array(self.control_points)
        
        ax.plot(points[:, 0], points[:, 1], points[:, 2], color='C0')
        #ax.scatter(control_points[:, 0], control_points[:, 1], control_points[:, 2])
        
        if show:
            fig.show()