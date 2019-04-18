import os
from logging import getLogger

import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
matplotlib.use('Agg')
import pylab as plt
from matplotlib.colors import LogNorm

from omnium import Analyser
from omnium.utils import get_cube, cm_to_inch

from scaffold.expt_settings import EXPT_DETAILS

logger = getLogger('scaf.precip_plot')


class PrecipPlot(Analyser):
    """Pick out precip timesteps and plot."""
    analysis_name = 'precip_plot'
    multi_expt = True
    input_dir = 'work/19700101T0000Z/{expt}_atmos'
    input_filename_glob = '{input_dir}/atmos.pp3.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.precip_plot.dummy']

    expts_to_plot = None

    def load(self):
        self.load_cubes()
        self.precips = {}
        for expt in self.task.expts:
            cubes = self.expt_cubes[expt]
            precip = get_cube(cubes, 4, 203)
            self.precips[expt] = precip
        output_dir = os.path.dirname(self.task.output_filenames[0])
        self.cached_hov_data = {}
        for expt in self.task.expts:
            cached_hov_x_fn = os.path.join(output_dir, '{}_cached_hov_x.nc'.format(expt))
            cached_hov_y_fn = os.path.join(output_dir, '{}_cached_hov_y.nc'.format(expt))
            if os.path.exists(cached_hov_x_fn) and os.path.exists((cached_hov_y_fn)):
                self.cached_hov_data[expt] = np.load(cached_hov_x_fn), np.load(cached_hov_y_fn)


    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        if not self.expts_to_plot:
            self.expts_to_plot = self.task.expts

        self._plot_snapshots()
        self._plot_hovmollers()

    def _plot_snapshots(self):
        # Need one of them to get shape.
        precip = list(self.precips.values())[0]

        max_precips = ['time_index, expt, max precip [mm/hr]']
        for i in range(precip.shape[0] - 100, precip.shape[0]):
            if len(self.expts_to_plot) == 10:
                fig = plt.figure(figsize=cm_to_inch(18, 10))
                gs = gridspec.GridSpec(3, 5,
                                       height_ratios=[1, 1, 0.6])
                axes = []
                for ax_index in range(len(self.expts_to_plot)):
                    if ax_index == 0:
                        ax = plt.subplot(gs[ax_index // 5, ax_index % 5])
                    else:
                        ax = plt.subplot(gs[ax_index // 5, ax_index % 5],
                                         sharex=axes[0], sharey=axes[0])
                    axes.append(ax)

                    if ax_index in [0, 1, 2, 3, 4]:
                        plt.setp(ax.get_xticklabels(), visible=False)
                    else:
                        ax.set_xlabel('x (km)')

                    if ax_index in [1, 2, 3, 4, 6, 7, 8, 9]:
                        plt.setp(ax.get_yticklabels(), visible=False)
                    else:
                        ax.set_ylabel('y (km)')

            else:
                gs = gridspec.GridSpec(2, len(self.expts_to_plot),
                                       height_ratios=[1, 0.2])
                axes = []
                for ax_index in range(len(self.expts_to_plot)):
                    fig = plt.figure(figsize=cm_to_inch(18, 8))
                    if ax_index == 0:
                        ax = plt.subplot(gs[0, ax_index])
                    else:
                        ax = plt.subplot(gs[0, ax_index], sharex=axes[0], sharey=axes[0])
                    axes.append(ax)
                    ax.set_xlabel('x (km)')

                    if ax_index != 0:
                        plt.setp(ax.get_yticklabels(), visible=False)
                    else:
                        ax.set_ylabel('y (km)')

            precip_max = 0
            for expt in self.expts_to_plot:
                precip = self.precips[expt]
                precip_max = max(precip[i].data.max(), precip_max)
                # vmax = precip_max * 3600
                vmax = 200

            for ax, expt in zip(axes, self.expts_to_plot):
                if expt in EXPT_DETAILS:
                    ucp_kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
                    ax.set_title(ucp_kwargs['label'])
                else:
                    ax.set_title(expt)

                precip = self.precips[expt]
                # precip in kg m-2 s-1, want mm hr-1:
                # /rho_water: m s-1
                # *1000: mm s-1
                # *3600: mm hr-1
                # N.B. rho_water = 1000 kg m-3.
                #import ipdb; ipdb.set_trace()
                precip_data = precip[i].data * 3600

                precip_min = 1e-4
                precip_vmin = 1e-1
                precip_data[precip_data < precip_min] = 0
                im = ax.imshow(precip_data, origin='lower',
                               interpolation='nearest', extent=[0, 256, 0, 256],
                               #vmin=0, vmax=precip_max * 3600)
                               norm=LogNorm(vmin=precip_vmin, vmax=vmax))

            for expt in self.task.expts:
                precip = self.precips[expt]
                precip_data = precip[i].data * 3600
                max_precips.append('{},{},{}'.format(i, expt, precip_data.max()))

            if len(self.expts_to_plot) == 10:
                plt.subplots_adjust(bottom=0.4)
                cbar_ax = fig.add_axes([0.1, 0.15, 0.8, 0.02])
            else:
                plt.subplots_adjust(bottom=0.3)
                cbar_ax = fig.add_axes([0.1, 0.2, 0.8, 0.04])
            plt.tight_layout()
            # cbar_ax = fig.add_axes([0.85, 0.27, 0.02, 0.46])
            cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
            # cbar.set_label('precip. (mm hr$^{-1}$)', rotation=270, labelpad=15)
            cbar.set_label('precip. (mm hr$^{-1}$)', labelpad=0)

            if i == 1849:
                # UCP figure!
                plt.savefig(self.file_path('UCP_time_index{}.png'.format(i)))
            plt.savefig(self.file_path('time_index{}.png'.format(i)))
            plt.close('all')
        self.save_text('max_precip.csv', '\n'.join(max_precips) + '\n')

    def _plot_hovmollers(self):
        if len(self.task.expts) != 4:
            logger.debug('only plot hovmollers for 4 expts')
            return

        # N.B. should divide into 20
        divider = 10

        fig = plt.figure(figsize=cm_to_inch(18, 15))
        # fig.subplots_adjust(bottom=0.15)

        # Displaying 2x4 data. Use 3x5 grid because need gap between left 2 and right 2 ax,
        # and need space for colorbar at bottom.
        gs = gridspec.GridSpec(3, 5,
                               height_ratios=[1, 1, 0.2],
                               width_ratios=[1, 1, 0.01, 1, 1])
        flattened_axes = []
        for i in [0, 1]:
            for j in [0, 1, 3, 4]:
                if i != 0 and j != 0:
                    flattened_axes.append(plt.subplot(gs[i, j], sharey=flattened_axes[0]))
                else:
                    flattened_axes.append(plt.subplot(gs[i, j]))
        colorbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.02])
        masked_hov = {}
        precip_min = 1e-1
        precip_max = 0
        for i, expt in enumerate(self.task.expts):
            precip = self.precips[expt]
            if expt in self.cached_hov_data:
                logger.debug('using cached data for {}', expt)
                masked_hov_x, masked_hov_y = self.cached_hov_data[expt]
            else:
                timesteps = precip.shape[0]
                hov_x = precip[timesteps - timesteps // divider:].data.mean(axis=1)
                hov_y = precip[timesteps - timesteps // divider:].data.mean(axis=2)

                # Mask out very small values of precip.
                mask_x = hov_x < 0.00001
                mask_y = hov_y < 0.00001

                masked_hov_x = np.ma.masked_array(hov_x, mask_x)
                masked_hov_y = np.ma.masked_array(hov_y, mask_y)

                output_dir = os.path.dirname(self.task.output_filenames[0])
                cached_hov_x_fn = os.path.join(output_dir, '{}_cached_hov_x.nc'.format(expt))
                cached_hov_y_fn = os.path.join(output_dir, '{}_cached_hov_y.nc'.format(expt))
                logger.debug('caching data for {}', expt)
                masked_hov_x.dump(cached_hov_x_fn)
                masked_hov_y.dump(cached_hov_y_fn)

            masked_hov[expt] = (masked_hov_x, masked_hov_y)
            precip_max = max(masked_hov_x.max(), masked_hov_y.max(), precip_max)

        for i, expt in enumerate(self.task.expts):
            (masked_hov_x, masked_hov_y) = masked_hov[expt]

            if expt in EXPT_DETAILS:
                expt_name, colour = EXPT_DETAILS[expt][0:2]
            else:
                expt_name = expt
                colour = 'k'

            ax1, ax2 = flattened_axes[i * 2], flattened_axes[i * 2 + 1]
            precip = self.precips[expt]

            start_time = precip.coord('time').points[0]
            times = precip.coord('time').points - start_time
            timelength_hours = times[-1] - times[0]
            timelength_days = timelength_hours / 24

            # Convert m to km.
            x = precip.coord('grid_longitude').points / 1e3
            y = precip.coord('grid_latitude').points / 1e3
            x_res = x[1] - x[0]
            y_res = y[1] - y[0]
            nx = len(x)
            ny = len(y)

            ax1.set_title('{} x-dir'.format(expt_name))
            im = ax1.imshow(masked_hov_x * 3600, interpolation='nearest',
                            extent=[0, x_res * nx, timelength_days, timelength_days - timelength_days / divider],
                            aspect='auto',
                            cmap=plt.get_cmap('Blues'),
                            #norm=LogNorm(), vmax=2e-3)
                            norm=LogNorm(vmin=precip_min, vmax=precip_max * 3600))
            ax1.set_xlabel('x (km)')

            ax2.set_title('{} y-dir'.format(expt_name))
            ax2.imshow(masked_hov_y * 3600, interpolation='nearest',
                       extent=[0, y_res * ny, timelength_days, timelength_days - timelength_days / divider],
                       aspect='auto',
                       cmap=plt.get_cmap('Blues'),
                       # norm=LogNorm(), vmax=2e-3)
                       norm=LogNorm(vmin=precip_min, vmax=precip_max * 3600))
            ax2.set_xlabel('y (km)')

            # Closure over expt
            def plot_phase_vel(x, y):
                dx = x[1] - x[0]
                dy = y[1] - y[0]
                # dx: km, dy: days
                phase_vel = dx / dy * 1e3 / 86400 # m/s
                logger.info('Phase velocity for {}: {} m/s', expt, phase_vel)
                ax1.plot(x, y, linestyle='--', color='k',
                         label='{:.1f} m s$^{{-1}}$'.format(phase_vel))
                ax1.legend(framealpha=1)
                return phase_vel

            x = [256, 0]
            if expt == 'S4W0Forced':
                y = [18.64, 19.5]
                plot_phase_vel(x, y)
            if expt == 'S4W5Forced':
                y = [19.125, 19.5]
                plot_phase_vel(x, y)
            if expt == 'S0W5Forced':
                y = [19.6, 19]
                plot_phase_vel(x, y)

            if i in [0, 2]:
                ax1.set_ylabel('time (days)')
            else:
                plt.setp(ax1.get_yticklabels(), visible=False)
            plt.setp(ax2.get_yticklabels(), visible=False)

        cb = plt.colorbar(im, cax=colorbar_ax, orientation='horizontal')
        cb.set_label('precip. (mm hr$^{-1}$)')

        plt.tight_layout()
        plt.show()
        plt.savefig(self.file_path('hovmollers.png'))

