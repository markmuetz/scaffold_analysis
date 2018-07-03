from logging import getLogger

import matplotlib

matplotlib.use('Agg')
import numpy as np
import iris

from cloud_tracking.utils import label_clds

from omnium.analyser import Analyser
from omnium.utils import get_cube
from omnium.consts import Re

from scaffold.scaffold_settings import settings

logger = getLogger('scaf.ca')


class CloudAnalyser(Analyser):
    """Performs thresholding of w, qcl to determine where clouds are.

    Does thresholding on height_levels.
    Performs analysis once for each combination of height_level, w_thresh, qcl_thresh.

    Outputs cloud_mask: result of combined (&'ed) w and qcl thresholding.
    cloud_mask has coords for the height_level, w_thresh and qcl_thresh used.
    Also output w, qcl and rho slices for future analysis.
    """
    analysis_name = 'cloud_analysis'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.pp1.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.{runid:03}.cloud_analysis.nc']
    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).pp1.nc'

    settings = settings

    def load(self):
        self.load_cubes()

    def run(self):
        self._apply_cloud_thresholds()
        self._label_clouds()

    def save(self, state, suite):
        self.save_results_cubes(state, suite)

    def _apply_cloud_thresholds(self):
        cubes = self.cubes

        w = get_cube(cubes, 0, 150)
        qcl = get_cube(cubes, 0, 254)
        rho = get_cube(cubes, 0, 253)

        rho.data = rho.data / Re**2
        rho.units = 'kg m-3'

        rho_height_levels = []
        for height_level in settings.height_levels:
            # N.B. this is right *if you are running using pp1*.
            rho_height_levels.extend([height_level, height_level + 1])

        # ONLY use these to get coords.
        w_slice = w[:, settings.height_levels]
        rho_slice = w[:, rho_height_levels]
        slice_coords = [(w.coord('time'), 0),
                        (w_slice.coord('model_level_number'), 1),
                        (w.coord('grid_latitude'), 2),
                        (w.coord('grid_longitude'), 3)]

        rho_coords = [(w.coord('time'), 0),
                      (rho_slice.coord('model_level_number'), 1),
                      (w.coord('grid_latitude'), 2),
                      (w.coord('grid_longitude'), 3)]

        # WARNING, you have to do it like this or it's VERY SLOW.
        # N.B. data deref'd *before* slice.
        logger.debug('slicing data')
        self.results['w_slice'] = iris.cube.Cube(w.data[:, settings.height_levels],
                                                 long_name='w_slice',
                                                 units='m s-1',
                                                 dim_coords_and_dims=slice_coords)
        theta_heights = w.coord('level_height').points[settings.height_levels]
        rho_heights = rho.coord('level_height').points[rho_height_levels]
        self.results['w_slice'].attributes['heights'] = theta_heights

        self.results['qcl_slice'] = iris.cube.Cube(qcl.data[:, settings.height_levels],
                                                   long_name='qcl_slice',
                                                   units='kg kg-1',
                                                   dim_coords_and_dims=slice_coords)
        self.results['qcl_slice'].attributes['heights'] = theta_heights

        self.results['rho_slice'] = iris.cube.Cube(rho.data[:, rho_height_levels],
                                                   long_name='rho_slice',
                                                   units='kg m-3',
                                                   dim_coords_and_dims=rho_coords)
        self.results['rho_slice'].attributes['heights'] = rho_heights

        # VERY SLOW if data you are reading are compressed.
        # w_slice = w[:, settings.height_levels]
        # self.results['w_slice'] = w_slice
        # self.results['w_slice'] = w_slice
        # self.results['qcl_slice'] = qcl[:, settings.height_levels].copy()
        # self.results['rho_slice'] = rho[:, rho_height_levels].copy()

        cloud_mask_data = np.zeros((w.shape[0],
                                    len(settings.height_levels),
                                    len(settings.w_threshs),
                                    len(settings.qcl_threshs),
                                    w.shape[-2],
                                    w.shape[-1]))

        for h_index, height_level in enumerate(settings.height_levels):
            for w_index, w_thresh in enumerate(settings.w_threshs):
                for qcl_index, qcl_thresh in enumerate(settings.qcl_threshs):
                    logger.debug('h, w, qcl index: ({}, {}, {})'.format(h_index, w_index, qcl_index))
                    w_mask = w[:, height_level].data > w_thresh
                    qcl_mask = qcl[:, height_level].data > qcl_thresh

                    cloud_mask_data[:, h_index, w_index, qcl_index] = (w_mask & qcl_mask).astype(int)

        w_thresh_coord = iris.coords.DimCoord(settings.w_threshs, long_name='w_thres', units='m s-1')
        qcl_thresh_coord = iris.coords.DimCoord(settings.qcl_threshs, long_name='qcl_thres', units='kg kg-1')
        cloud_mask_cube = iris.cube.Cube(cloud_mask_data,
                                         long_name='cloud_mask',
                                         dim_coords_and_dims=[(w.coord('time'), 0),
                                                              (w_slice.coord('model_level_number'), 1),
                                                              (w_thresh_coord, 2),
                                                              (qcl_thresh_coord, 3),
                                                              (w.coord('grid_latitude'), 4),
                                                              (w.coord('grid_longitude'), 5)])

        self.results['cloud_mask'] = cloud_mask_cube

    def _label_clouds(self):
        cloud_mask_cube = self.results['cloud_mask']
        w_slice = self.results['w_slice']

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        qcl_thresh_coord = cloud_mask_cube.coord('qcl_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')

        # level_number refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, level_number in enumerate(level_number_coord.points):

            for thresh_index in range(w_thresh_coord.shape[0]):
                # N.B. I just take the diagonal indices.
                w_thresh = w_thresh_coord.points[thresh_index]
                qcl_thresh = qcl_thresh_coord.points[thresh_index]

                labelled_clouds_cube = w_slice[:, height_level_index].copy()
                labelled_clouds_cube.units = ''

                labelled_clouds_data = np.zeros_like(labelled_clouds_cube.data)
                for time_index in range(cloud_mask_cube.data.shape[0]):
                    cloud_mask_ss = cloud_mask_cube[time_index,
                                                    height_level_index,
                                                    thresh_index,
                                                    thresh_index].data.astype(bool)
                    max_cld_index, labelled_clouds = label_clds(cloud_mask_ss, diagonal=True)
                    labelled_clouds_data[time_index] = labelled_clouds

                labelled_clouds_cube_id = 'labelled_clouds_z{}_w{}_qcl{}'.format(level_number,
                                                                                 w_thresh,
                                                                                 qcl_thresh)

                labelled_clouds_cube.rename(labelled_clouds_cube_id)
                labelled_clouds_cube.data = labelled_clouds_data
                self.results[labelled_clouds_cube_id] = labelled_clouds_cube
