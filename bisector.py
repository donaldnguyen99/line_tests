# Method
# 1. Use left side points and spline interpolation to get the right side points
# 2. Use right side points and spline interpolation to get the left side points
# 3. Use left side points and right side points to get the midpoints

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, minimize


class Bisector:

    def __init__(self, x=None, y=None):
        # if no x or y, create a fake absorption spectrum
        if x is None:
            self.x = np.linspace(0, 10, 15)
        else:
            self.x = x
        if y is None:
            self.y = 1 - self.gaussian(self.x, 5.5, 1) - self.gaussian(self.x, 5.7, 0.2) * 0.2
        else:
            self.y = y
        self.spline = interp1d(self.x, self.y, kind='cubic')
        self.id_min = np.argmin(self.y)
        self.interp_min_x = minimize(self.spline, self.x[self.id_min], method='nelder-mead').x[0]
        self.interp_min_y = self.spline(self.interp_min_x)
        self.bis_max_y = 0.9
        self.side_points = self.get_points_on_side()
        self.line_bisectors = np.concatenate((self.midpoints_from('left')[0], self.midpoints_from('right')[0]), axis=0)
        self.sorted_line_bisectors = self.line_bisectors[self.line_bisectors[:, 1].argsort()] # sort by y
        self.spline_bisectors = interp1d(self.sorted_line_bisectors[:, 1], self.sorted_line_bisectors[:, 0], kind='cubic') # interp func for bisector's x given y


    @staticmethod
    @np.vectorize
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


    def solve_for_x(self, y, x_guess=0, func=None):
        if func is None:
            func = self.spline    
        return fsolve(lambda x: func(x) - y, x_guess)
    

    def get_points_on_side(self, y_max=None):
        if y_max is None:
            y_max = self.bis_max_y
        leftside_max_ind = self.id_min
        rightside_min_ind = self.id_min
        if self.x[self.id_min] <= self.interp_min_x:
            leftside_max_ind += 1
            rightside_min_ind += 1
        elif self.x[self.id_min] == self.interp_min_x:
            rightside_min_ind += 1
        left = np.vstack((self.x[:leftside_max_ind], self.y[:leftside_max_ind])).T
        right = np.vstack((self.x[rightside_min_ind:], self.y[rightside_min_ind:])).T
        points_dict = {
            'left': left.compress(left[:, 1] <= y_max, axis=0),
            'right': right.compress(right[:, 1] <= y_max, axis=0)
        }
        return points_dict
            

    def midpoints_from(self, side):
        interp_side = np.empty(self.side_points[side].shape)
        midpoints = np.empty(self.side_points[side].shape)
        for i, coord in enumerate(self.side_points[side]):
            interp_side[i][0] = self.solve_for_x(coord[1], self.interp_min_x + (self.interp_min_x - coord[0]), self.spline)
            interp_side[i][1] = coord[1]
            midpoints[i][0] = (coord[0] + interp_side[i][0]) / 2
            midpoints[i][1] = coord[1]
        return midpoints, interp_side
            


    def plot(self):
        # plot fake absorption spectra
        plt.figure()

        new_x = np.linspace(0, 10, 100)
        plt.plot(new_x, self.spline(new_x), label='cubic spline')
        plt.scatter(self.x, self.y, label='data')

        # plot side points and corresponding bisectors
        plt.scatter(self.side_points['left'][:, 0], self.side_points['left'][:, 1], marker='s', c='orange', label='left side points')
        plt.scatter(self.midpoints_from('left')[0][:, 0], self.midpoints_from('left')[0][:, 1], marker='*', c='orange', label='left midpoints')
        plt.scatter(self.side_points['right'][:, 0], self.side_points['right'][:, 1], marker='^', c='green', label='right side points')
        plt.scatter(self.midpoints_from('right')[0][:, 0], self.midpoints_from('right')[0][:, 1], marker='*', c='green', label='right midpoints')

        # plot interpolation function for bisectors
        plt.plot(self.spline_bisectors(np.linspace(self.sorted_line_bisectors[0, 1], self.sorted_line_bisectors[-1, 1], 100)), np.linspace(self.sorted_line_bisectors[0, 1], self.sorted_line_bisectors[-1, 1], 100), label='interpolation function for bisectors')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    bisector = Bisector()
    bisector.plot()