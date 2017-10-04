import os
from logging import getLogger

import numpy as np

import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr
from cloud_tracking.utils import label_clds
from cloud_tracking import Tracker
from cloud_tracking.cloud_tracking_analysis import generate_stats, plot_stats, output_stats_to_file

logger = getLogger('scaf.cta')


class CloudTrackAnalyser(Analyser):
    analysis_name = 'cloud_track_analysis'
    multi_file = True

    def load(self):
        self.append_log('Override load')

        cloud_mask_cubes = []
        cloud_mask_id = 'cloud_mask'
        for filename in self.filenames:
            cubes = iris.load(filename)
            cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)
            cloud_mask_cubes.append(cloud_mask_cube)

        self.cubes = iris.cube.CubeList(cloud_mask_cubes).concatenate()
        self.append_log('Override loaded')

    def run_analysis(self):
        cubes = self.cubes

        self.trackers = {}
        self.all_stats = {}
        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')
        logger.debug(cloud_mask_cube.shape)

        # height_level refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, height_level in enumerate(level_number_coord.points):
            for thresh_index in range(w_thresh_coord.shape[0]):
                cld_field = np.zeros(cloud_mask_cube[:, height_level_index, thresh_index, thresh_index].shape, dtype=int)
                cld_field_cube = cloud_mask_cube[:, height_level_index, thresh_index, thresh_index].copy()
                cld_field_cube.rename('cloud_field')

                for time_index in range(cloud_mask_cube.shape[0]):
                    cloud_mask_ss = cloud_mask_cube[time_index,
                                                    height_level_index,
                                                    thresh_index,
                                                    thresh_index].data.astype(bool)
                    max_label, cld_labels = label_clds(cloud_mask_ss, True)
                    cld_field[time_index] = cld_labels
                cld_field_cube.data = cld_field

                tracker = Tracker(cld_field_cube.slices_over('time'), store_working=True)
                tracker.track()
                tracker.group()

                stats = generate_stats(self.expt, tracker)
                self.trackers[(height_level_index, thresh_index)] = tracker
                self.all_stats[(height_level_index, thresh_index)] = stats

    def display_results(self):
        self.append_log('displaying results')
        figpath = self.figpath('cloud_tracking')

        cubes = self.cubes
        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)
        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')
        for height_level_index, height_level in enumerate(level_number_coord.points):
            for thresh_index in range(w_thresh_coord.shape[0]):
                stats = self.all_stats[(height_level_index, thresh_index)]
                tracker = self.trackers[(height_level_index, thresh_index)]

                filename = 'atmos.cloud_tracking_z{}_t{}.'.format(height_level_index, thresh_index)

                plot_stats(self.expt, os.path.dirname(figpath), filename, [stats])
                output_stats_to_file(self.expt, os.path.dirname(figpath), filename + 'txt', tracker, stats)
