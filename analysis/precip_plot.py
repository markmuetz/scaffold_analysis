import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
import iris
from matplotlib.colors import LogNorm

from omnium.analyzer import Analyzer
from omnium.utils import get_cube
from omnium.consts import L


class PrecipPlot(Analyzer):
    """Pick out precip timesteps and plot."""
    analysis_name = 'precip_plot'
    multi_expt = True
    expts_to_plot = ['S0', 'S4']

    def _plot(self):
        precips = {}
        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            precip = get_cube(cubes, 4, 203)
            precips[expt] = precip

        max_precips = ['time_index, expt, max precip [mm/hr]']
        for i in range(precip.shape[0] - 100, precip.shape[0]):
            fig, axes = plt.subplots(1, len(self.expts_to_plot))

            precip_max = 0
            for expt in self.expts_to_plot:
                precip = precips[expt]
                precip_max = max(precip[i].data.max(), precip_max)

            for ax, expt in zip(axes, self.expts_to_plot):
                ax.set_title(expt)
                if expt == self.expts[0]:
                    ax.set_ylabel('y (km)')
                else:
                    ax.get_yaxis().set_visible(False)
                ax.set_xlabel('x (km)')

                precip = precips[expt]
                # precip in kg m-2 s-1, want mm hr-1:
                # /rho_water: m s-1
                # *1000: mm s-1
                # *3600: mm hr-1
                # N.B. rho_water = 1000 kg m-3.
                #import ipdb; ipdb.set_trace()
                precip_data = precip[i].data * 3600

                precip_min = 1e-4
                precip_data[precip_data < precip_min] = 0
                im = ax.imshow(precip_data, origin='lower', 
                               interpolation='nearest', extent=[0, 256, 0, 256],
                               #vmin=0, vmax=precip_max * 3600)
                               norm=LogNorm(vmin=precip_min, vmax=precip_max * 3600))

            for expt in self.expts:
                precip = precips[expt]
                precip_data = precip[i].data * 3600
                max_precips.append('{},{},{}'.format(i, expt, precip_data.max()))

            plt.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.89, 0.27, 0.02, 0.46])
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label('rainfall (mm hr$^{-1}$)', rotation=270, labelpad=20)

            plt.savefig(self.figpath('time_index{}.png'.format(i)))
            plt.close('all')
        self.save_text('max_precip.csv', '\n'.join(max_precips) + '\n')

    def run_analysis(self):
        pass

    def display_results(self):
        """Save all results for surf flux analysis."""
        self._plot()

