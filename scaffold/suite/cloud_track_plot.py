import os
from logging import getLogger
import pickle

import numpy as np
import matplotlib.pyplot as plt

from scaffold.expt_settings import EXPT_DETAILS

from omnium import Analyser
from omnium.utils import cm_to_inch

logger = getLogger('scaf.ctp')

def plot_pdf_for_expt(ax, expt_name, all_stats):
    plt_names = ['all_lifetimes', 'linear_lifetimes', 'nonlinear_lifetimes']


    combined_plot_kwargs = {
        'all_lifetimes': {'label': 'all ', 'color': 'grey'},
        'linear_lifetimes': {'label': 'simple', 'color': 'g'},
        'nonlinear_lifetimes': {'label': 'complex', 'color': 'r'},
    }

    for stats in all_stats:
        # Note - normalization done w.r.t. all_lifetimes so that all figs have equiv axes.
        all_lifetimes_hist_sum = np.histogram(stats['all_lifetimes'],
                                              bins=80, range=(0, 400))[0].sum()
        max_height = 0
        for plt_name in plt_names:
            lifetimes = stats[plt_name]
            if not len(lifetimes):
                logger.warning('No lifetimes for {}'.format(plt_name))
                continue

            hist, bins = np.histogram(lifetimes, bins=80, range=(0, 400))
            widths = bins[1:] - bins[:-1]
            centres = (bins[1:] + bins[:-1]) / 2
            heights = (hist / all_lifetimes_hist_sum) / widths
            max_height = max(max_height, heights.max())

            ax.plot(centres, heights * widths, **combined_plot_kwargs[plt_name])

    ax.set_xlim((0, 400))
    ax.set_ylim((0, 0.1))


class CloudTrackPlot(Analyser):
    """Tracks clouds using method similar to Plant 2009.

    Change is due to temporal resolution of data being less, take account of this by first
    calculating spatial correlation then using this to project the cloud field at a given height
    forward in time. Most heavy lifting is handled by cloud_tracking package."""
    analysis_name = 'cloud_track_plot'
    multi_expt = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename = '{input_dir}/atmos.cloud_track_analysis.all_stats.pkl'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.cloud_track_plot.dummy']

    def load(self):
        self.all_stats = {}
        for expt, fn in zip(self.task.expts, self.task.filenames):
            self.append_log('loading {} for {}'.format(fn, expt))
            with open(fn, 'rb') as f:
                self.all_stats[expt] = pickle.load(f)

    def run(self):
        self.thresh = 2  # 0, 1 or 2 -- 10% lower, actual, 10% higher.
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        if not len(self.task.expts) == 4:
            logger.debug('can only run with 4 expts')
            return

        fig, axes = plt.subplots(2, 2, figsize=cm_to_inch(16, 14))

        height_level_index, thresh_index = (1, self.thresh)
        for i, (expt, ax) in enumerate(zip(self.task.expts, axes.flatten())):
            all_stats = self.all_stats[expt]

            stats = all_stats[(1, self.thresh)]

            if expt in EXPT_DETAILS:
                expt_name = EXPT_DETAILS[expt][0]
            else:
                expt_name = expt

            plot_pdf_for_expt(ax, expt_name, [stats])

            if i in [0, 2]:
                ax.set_ylabel('Frequency of lifecycle')
            else:
                plt.setp(ax.get_yticklabels(), visible=False)

            if i in [2, 3]:
                ax.set_xlabel('Lifetime (min)')
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            ax.set_title('{}'.format(expt_name))
            if i == 1:
                ax.legend(loc='upper right')


        plt.tight_layout()
        plt.savefig(self.file_path('cloud_tracking_z{}_t{}'.format(height_level_index,
                                                                   thresh_index)))
