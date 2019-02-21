import os
from logging import getLogger
import pickle

import numpy as np
from cloud_tracking import Tracker
from cloud_tracking.cloud_tracking_analysis import generate_stats, plot_stats, output_stats_to_file

from omnium import Analyser
from omnium.utils import get_cube_from_attr

logger = getLogger('scaf.ctp')


class CloudTrackPlot(Analyser):
    """Tracks clouds using method similar to Plant 2009.

    Change is due to temporal resolution of data being less, take account of this by first
    calculating spatial correlation then using this to project the cloud field at a given height
    forward in time. Most heavy lifting is handled by cloud_tracking package."""
    analysis_name = 'cloud_track_plot'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filenames = ['{input_dir}/atmos.cloud_track_analysis.all_stats.pkl',
                       '{input_dir}/atmos.cloud_track_analysis.trackers.pkl', ]
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.cloud_track_plot.dummy']

    def load(self):
        with open(self.task.filenames[0], 'rb') as f:
            self.all_stats = pickle.load(f)
        with open(self.task.filenames[1], 'rb') as f:
            self.trackers = pickle.load(f)

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self.append_log('displaying results')
        figpath = self.file_path('cloud_tracking')

        for tracker_key in self.trackers.keys():
            height_level_index, thresh_index = tracker_key
            stats = self.all_stats[tracker_key]
            tracker = self.trackers[tracker_key]
            filename = 'atmos.cloud_tracking_z{}_t{}.'.format(height_level_index, thresh_index)

            plot_stats(self.task.expt, os.path.dirname(figpath), filename, [stats])
            output_stats_to_file(self.task.expt, os.path.dirname(figpath), filename + 'txt',
                                 tracker, stats)
