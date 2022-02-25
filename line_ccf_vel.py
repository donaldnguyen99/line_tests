import sys
import path
import os
sys.path.append(os.path.abspath(__file__)+'../expres_pipeline')
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
import matplotlib.pyplot as plt

class CCF_Vel_Calibration:

    def __init__(self, x=None, y=None, x_template=None, y_template=None):
        # Initialize x and y
        if x is not None:
            if y is None:
                raise ValueError('y must be provided if x is provided')
            self.x = x
        
        if y is not None:
            if x is None:
                raise ValueError('x must be provided if y is provided')
            self.y = y
        
        if x is None and y is None:
            self.x, self.y = self.example_line_profile(0.2, 1, limits=5, num_points=100)

        # Initialize x_template and y_template
        if x_template is not None:
            if y_template is None:
                raise ValueError('y_template must be provided if x_template is provided')
            self.x_template = x_template
        
        if y_template is not None:
            if x_template is None:
                raise ValueError('x_template must be provided if y_template is provided')
            self.y_template = y_template
        
        if x_template is None and y_template is None:
            self.x_template, self.y_template = self.example_line_profile(limits=10, num_points=200)
        
        
        # self.x_interp = np.linspace(x.min(), x.max(), 1000)

    def gaussian(self, x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    def example_line_profile(self, mu=0, sig=1, limits=5, num_points=10):
        x = np.linspace(-limits * sig + mu, limits * sig + mu, num_points)
        y =  1 - self.gaussian(x, mu, sig)
        return x, y

    def get_ccf(self):
        pass

    def filter_same_range(self, x, other_x):
        return (x >= other_x.min()) & (x <= other_x.max())

    def plot(self):
        # x, y, x_template, and y_template are initialized in __init__

        # Line is used to generate the interpolation function
        line_interp = interp1d(self.x, self.y, kind='cubic')
        template_interp = interp1d(self.x_template, self.y_template, kind='cubic')

        # need to filter out points that are outside the range of the line
        #    ensure x_template has a large enough range to cover the range of x
        # assert(self.x_template.max() >= self.x.max() and self.x_template.min() <= self.x.min())

        line_interp_x = np.linspace(self.x.min(), self.x.max(), len(self.x) * 1000)
        line_interp_y = line_interp(line_interp_x)
        template_interp_x = np.linspace(self.x_template.min(), self.x_template.max(), len(self.x_template) * 1000)
        template_interp_y = template_interp(template_interp_x)
        # same_range_filter = (self.x_template >= self.x.min()) & (self.x_template <= self.x.max())
        # line_x_interp = self.x_template[same_range_filter]
        # line_y_interp = line_interp(line_x_interp)

        ccf_input1 = template_interp_y
        ccf_input2 = line_interp_y
        ccf = correlate(ccf_input1, ccf_input2, "full")
        ccf_x_ind = correlation_lags(len(ccf_input1), len(ccf_input2), "full")
        lag = ccf_x_ind[np.argmax(ccf)]
        wave_shift = line_interp_x[lag]

        # Plotting
        fig, axs = plt.subplots(3, 1, sharex=False)
        axs[0].plot(self.x_template, self.y_template, label='template')
        axs[0].plot(self.x, self.y, label='line')
        axs[0].plot(self.x - wave_shift, self.y, label='line shifted')
        axs[0].legend()
        
        # multiply by a factor that shifts the center (ratio of old to new)
        # multiply wave by z + 1
        axs[1].plot(template_interp_x, template_interp_y, label='template_interp')
        # axs[1].plot(line_interp_x, line_interp_y, label='line_interp')
        axs[1].plot(line_interp_x - wave_shift, line_interp_y, label='line_interp shifted')
        axs[1].scatter(template_interp_x[np.argmin(template_interp_y)], template_interp_y[np.argmin(template_interp_y)])
        # axs[1].scatter(line_interp_x[np.argmin(line_interp_y)], line_interp_y[np.argmin(line_interp_y)])
        axs[1].scatter((line_interp_x - wave_shift)[np.argmin(line_interp_y)], line_interp_y[np.argmin(line_interp_y)])
        axs[1].legend()
        # axs[0].scatter(line_x_interp, line_y_interp, label='interp_line')

        # axs[2].plot(np.arange(len(ccf)), ccf, label='ccf')
        axs[2].plot(ccf_x_ind, ccf, label='ccf_lags')
        axs[2].legend()
        plt.show(block=True)

if __name__ == "__main__":
    ccf = CCF_Vel_Calibration()
    ccf.plot()