import os
from logging import getLogger

import numpy as np
from cloud_tracking import Tracker
from cloud_tracking.cloud_tracking_analysis import generate_stats, plot_stats, output_stats_to_file

from omnium import Analyser
from omnium.utils import get_cube_from_attr

logger = getLogger('scaf.cta')


class CloudTrackAnalyser(Analyser):
    """Tracks clouds using method similar to Plant 2009.

    Change is due to temporal resolution of data being less, take account of this by first
    calculating spatial correlation then using this to project the cloud field at a given height
    forward in time. Most heavy lifting is handled by cloud_tracking package."""
    analysis_name = 'cloud_track_analysis'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.cloud_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.cloud_track_analysis.dummy']
    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).cloud_analysis.nc'
    min_runid = 480
    # No max.
    # max_runid = 308

    def load(self):
        self.load_cubes()

    def run(self):
        cubes = self.cubes

        self.trackers = {}
        self.all_stats = {}
        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        qcl_thresh_coord = cloud_mask_cube.coord('qcl_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')
        logger.debug(cloud_mask_cube.shape)

        # height_level refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, level_number in enumerate(level_number_coord.points):
            for thresh_index in range(w_thresh_coord.shape[0]):

                logger.debug('height_index, thresh_index: {}, {}'.format(height_level_index,
                                                                         thresh_index))

                w_thresh = w_thresh_coord.points[thresh_index]
                qcl_thresh = qcl_thresh_coord.points[thresh_index]
                labelled_clouds_cube_id = 'labelled_clouds_z{}_w{}_qcl{}'.format(level_number,
                                                                                 w_thresh,
                                                                                 qcl_thresh)
                labelled_clouds_cube = get_cube_from_attr(cubes,
                                                          'omnium_cube_id',
                                                          labelled_clouds_cube_id)
                cld_field = np.zeros(cloud_mask_cube[:, height_level_index, thresh_index, thresh_index].shape, dtype=int)
                cld_field_cube = cloud_mask_cube[:, height_level_index, thresh_index, thresh_index].copy()
                cld_field_cube.rename('cloud_field')

                for time_index in range(cloud_mask_cube.shape[0]):
                    labelled_clouds_ss = labelled_clouds_cube[time_index].data.astype(int)
                    cld_field[time_index] = labelled_clouds_ss
                cld_field_cube.data = cld_field

                tracker = Tracker(cld_field_cube.slices_over('time'), store_working=True)
                tracker.track()
                tracker.group()

                stats = generate_stats(self.task.expt, tracker)
                self.trackers[(height_level_index, thresh_index)] = tracker
                self.all_stats[(height_level_index, thresh_index)] = stats

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self.append_log('displaying results')
        figpath = self.file_path('cloud_tracking')

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

                plot_stats(self.task.expt, os.path.dirname(figpath), filename, [stats])
                output_stats_to_file(self.task.expt, os.path.dirname(figpath), filename + 'txt', tracker, stats)
