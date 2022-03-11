import sys
import path
import os
sys.path.append(os.path.abspath(__file__)+'../expres_pipeline/expres')
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, minimize
from scipy.signal import correlate, correlation_lags
import matplotlib.pyplot as plt

# import expres.rv.ccf.fit_rv as fit_rv

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

    @staticmethod
    def gaussian(x, mu, sig, A=1):
        return A * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    def example_line_profile(self, mu=0, sig=1, limits=5, num_points=10):
        x = np.linspace(-limits * sig + mu, limits * sig + mu, num_points)
        y =  1 - self.gaussian(x, mu, sig)
        return x, y

    def get_wave_shift(self):
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
        template_core_x = minimize(template_interp, template_interp_x[np.argmin(template_interp_y)]).x[0]
        # same_range_filter = (self.x_template >= self.x.min()) & (self.x_template <= self.x.max())
        # line_x_interp = self.x_template[same_range_filter]
        # line_y_interp = line_interp(line_x_interp)

        ccf_input1 = template_interp_y
        ccf_input2 = line_interp_y
        # smaller range instead of full
        ccf = correlate(ccf_input1, ccf_input2, "valid")
        ccf_x_ind = correlation_lags(len(ccf_input1), len(ccf_input2), "valid")
        lag_index = np.argmax(ccf)
        lag = ccf_x_ind[lag_index]

        # fit gaussian to ccf centered at lag
        popt_full, _  = curve_fit(self.gaussian, ccf_x_ind, ccf, p0=[lag, ccf_x_ind[-1] - lag, ccf[lag_index]])
        # ccf_gen_fit = self.gaussian(ccf_x_ind, *popt_full)

        best_sig = popt_full[1] / 10
        add_side_points = 200
        ccf_seg_x = ccf_x_ind[lag_index-add_side_points:lag_index+add_side_points+1]
        ccf_seg_y = ccf[lag_index-add_side_points:lag_index+add_side_points+1]
        popt, pcov = curve_fit(self.gaussian, ccf_seg_x, ccf_seg_y, p0=[popt_full[0], best_sig, popt_full[2]])
        mu, sig, amp = popt[0], popt[1], popt[2]
        # print('mu = {}, sig = {}, amp = {}'.format(mu, sig, amp))
        # ccf_fitted_seg = self.gaussian(ccf_seg_x, *popt)
        # ccf_fitted = self.gaussian(ccf_x_ind, *popt)

        if mu < lag:
            start_index = lag_index
            while (mu < ccf_x_ind[start_index]):
                start_index -= 1
                end_index = start_index + 1
        elif mu > lag:
            end_index = lag_index
            while (mu > ccf_x_ind[end_index]):
                end_index += 1
            start_index = end_index - 1

        # TODO - get wave shift error from ccf sig
        
        interp_perc_between_indices = (mu - ccf_x_ind[start_index]) / (ccf_x_ind[end_index] - ccf_x_ind[start_index])
        wave_shift = interp_perc_between_indices * (line_interp_x[ccf_x_ind[end_index]] - line_interp_x[ccf_x_ind[start_index]]) + line_interp_x[ccf_x_ind[start_index]]
        return wave_shift, template_core_x

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

        line_interp_x = np.linspace(self.x.min(), self.x.max(), len(self.x) * 100)
        line_interp_y = line_interp(line_interp_x)
        template_interp_x = np.linspace(self.x_template.min(), self.x_template.max(), len(self.x_template) * 100)
        template_interp_y = template_interp(template_interp_x)
        # same_range_filter = (self.x_template >= self.x.min()) & (self.x_template <= self.x.max())
        # line_x_interp = self.x_template[same_range_filter]
        # line_y_interp = line_interp(line_x_interp)

        ccf_input1 = template_interp_y
        ccf_input2 = line_interp_y
        # smaller range instead of full
        ccf = correlate(ccf_input1, ccf_input2, "valid")
        ccf_x_ind = correlation_lags(len(ccf_input1), len(ccf_input2), "valid")
        lag_index = np.argmax(ccf)
        lag = ccf_x_ind[lag_index]

        # fit gaussian to ccf centered at lag
        popt_full, _  = curve_fit(self.gaussian, ccf_x_ind, ccf, p0=[lag, ccf_x_ind[-1] - lag, 65000])
        ccf_gen_fit = self.gaussian(ccf_x_ind, *popt_full)

        best_sig = popt_full[1] / 10
        add_side_points = 20
        ccf_seg_x = ccf_x_ind[lag_index-add_side_points:lag_index+add_side_points+1]
        ccf_seg_y = ccf[lag_index-add_side_points:lag_index+add_side_points+1]
        my_gaussian = lambda x, mu, sig: popt_full[2] * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
        popt, pcov = curve_fit(self.gaussian, ccf_seg_x, ccf_seg_y, p0=[popt_full[0], best_sig, popt_full[2]])
        mu, sig, amp = popt[0], popt[1], popt[2]
        print('mu = {}, sig = {}'.format(mu, sig))
        ccf_fitted_seg = self.gaussian(ccf_seg_x, *popt)
        ccf_fitted = self.gaussian(ccf_x_ind, *popt)

        wave_shift = mu - lag
        if mu < lag:
            start_index = lag_index
            while (mu < ccf_x_ind[start_index]):
                start_index -= 1
                end_index = start_index + 1
        elif mu > lag:
            end_index = lag_index
            while (mu > ccf_x_ind[end_index]):
                end_index += 1
            start_index = end_index - 1
        
        interp_perc_between_indices = (mu - ccf_x_ind[start_index]) / (ccf_x_ind[end_index] - ccf_x_ind[start_index])
        wave_shift = interp_perc_between_indices * (line_interp_x[ccf_x_ind[end_index]] - line_interp_x[ccf_x_ind[start_index]]) + line_interp_x[ccf_x_ind[start_index]]
        # wave_shift = line_interp_x[lag]

        # Plotting
        fig, axs = plt.subplots(4, 1, sharex=False)
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
        axs[2].plot(ccf_x_ind, ccf_gen_fit, label='ccf general fit')
        axs[2].plot(ccf_x_ind, ccf_fitted, label='ccf fit')
        axs[2].legend()

        axs[3].plot(ccf_seg_x, ccf_seg_y, label='ccf original')
        axs[3].plot(ccf_seg_x, ccf_fitted_seg, label='ccf_fitted')
        axs[3].legend()
        plt.show(block=True)


        # ccf fit gaussian at top of the ccf peak to find and interp to wavelength
        # v_grid = ccf
        # rv = fit_rv(v_grid, ccf, np.zeroes(ccf.shape))
        # print(rv)

if __name__ == "__main__":
    ccf = CCF_Vel_Calibration()
    ccf.plot()