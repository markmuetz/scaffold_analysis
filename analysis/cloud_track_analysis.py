import numpy as np

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr
from cloud_tracking.utils import label_clds
from cloud_tracking import Tracker


class CloudTrackAnalyzer(Analyzer):
    analysis_name = 'cloud_track_analysis'

    def run_analysis(self):
        cubes = self.cubes

        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')

        # height_level refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, height_level in enumerate(level_number_coord.points):
            for thresh_index in range(w_thresh_coord.shape[0]):
                cld_field = np.zeros(cloud_mask_cube.shape, dtype=int)
                cld_field_cube = cloud_mask_cube[:, height_level_index, thresh_index, thresh_index].copy()
                cld_field_cube.rename('cloud_field')

                for time_index in range(cloud_mask_cube.data.shape[0]):

                    cloud_mask_ss = cloud_mask_cube[time_index,
                                                    height_level_index,
                                                    thresh_index,
                                                    thresh_index].data.astype(bool)
                    max_label, cld_labels = label_clds(cloud_mask_ss, True)
                    cld_field[time_index] = cld_labels
                cld_field_cube.data = cld_field

                tracker = Tracker(cld_field_cube.slices_over('time'))
                tracker.track()
                # proj_cld_field_cube = cld_field_cube.copy()
                # proj_cld_field_cube.data = tracker.proj_cld_field.astype(float)
                # iris.save(proj_cld_field_cube, 'output/{}_proj_cld_field.nc'.format(expt))

                tracker.group()
