import sys
import path
sys.path.append(path.abspath(__file__)+'../expres_pipeline')
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class CCF_Vel_Calibration:

    def __init__(self, x=np.empty(0), y=np.empty(0), x_template=None, y_template=None):
        self.x = x
        self.y = y
        self.x_template = x_template
        self.y_template = y_template
        # self.x_interp = np.linspace(x.min(), x.max(), 1000)

    def gaussian(self, x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    def example_line_profile(self, mu=0, sig=1):
        x = np.linspace(-5 * sig + mu, 5 * sig + mu, 1000)
        y =  1 - self.gaussian(x, mu, sig)
        return x, y

    def get_ccf(self):
        pass

    def plot(self):
        
        line_template_x, line_template_y = self.example_line_profile()
        line_x, line_y = self.example_line_profile(0.3, 0.8)
        line_interp = interp1d(line_x, line_y, kind='cubic')
        # need to filter out points that are outside the range of the line
        #    ensure line_x has a large enough range
        # multiply wave by z + 1
        line_y_interp = line_interp(line_template_x)
        # ccf = np.correlate(line_y, line_y_interp, "full")
        plt.figure()
        plt.plot(line_template_x, line_template_y)
        plt.plot(line_x, line_y)
        plt.plot(line_template_x, line_y_interp)
        plt.show()

if __name__ == "__main__":
    ccf = CCF_Vel_Calibration()
    ccf.plot()