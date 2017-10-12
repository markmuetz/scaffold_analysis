from logging import getLogger

import matplotlib

matplotlib.use('Agg')
import numpy as np
import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube
from omnium.consts import Re
from omnium.omnium_errors import OmniumError

logger = getLogger('scaf.ca')


class CloudAnalyser(Analyser):
    analysis_name = 'cloud_analysis'
    single_file = True

    def set_config(self, config):
        super(CloudAnalyser, self).set_config(config)
        self.height_levels = [int(l) for l in config['height_levels'].split(',')]
        self.w_threshs = [float(t) for t in config['w_threshs'].split(',')]
        self.qcl_threshs = [float(t) for t in config['qcl_threshs'].split(',')]
        if len(self.w_threshs) != len(self.qcl_threshs):
            raise OmniumError('w_threshs and qcl_threshs must be same length')
            
    def run_analysis(self):
        cubes = self.cubes

        w = get_cube(cubes, 0, 150)
        qcl = get_cube(cubes, 0, 254)
        rho = get_cube(cubes, 0, 253)

        rho.data = rho.data / Re**2
        rho.units = 'kg m-3'

        rho_height_levels = []
        for height_level in self.height_levels:
            # N.B. this is right *if you are running using pp1*.
            rho_height_levels.extend([height_level, height_level + 1])

        # ONLY use these to get coords.
        w_slice = w[:, self.height_levels]
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
        self.results['w_slice'] = iris.cube.Cube(w.data[:, self.height_levels],
                                                 long_name='w_slice',
                                                 units='m s-1',
                                                 dim_coords_and_dims=slice_coords)
        theta_heights = w.coord('level_height').points[self.height_levels]
        rho_heights = rho.coord('level_height').points[rho_height_levels]
        self.results['w_slice'].attributes['heights'] = theta_heights
                                                 
        self.results['qcl_slice'] = iris.cube.Cube(qcl.data[:, self.height_levels],
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
        # w_slice = w[:, self.height_levels]
        # self.results['w_slice'] = w_slice
        # self.results['w_slice'] = w_slice
        # self.results['qcl_slice'] = qcl[:, self.height_levels].copy()
        # self.results['rho_slice'] = rho[:, rho_height_levels].copy()

        cloud_mask_data = np.zeros((w.shape[0], 
                                    len(self.height_levels),
                                    len(self.w_threshs),
                                    len(self.qcl_threshs),
                                    w.shape[-2],
                                    w.shape[-1]))
                              
        for h_index, height_level in enumerate(self.height_levels):
            for w_index, w_thresh in enumerate(self.w_threshs):
                for qcl_index, qcl_thresh in enumerate(self.qcl_threshs):
                    logger.debug('h, w, qcl index: ({}, {}, {})'.format(h_index, w_index, qcl_index))
                    w_mask = w[:, height_level].data > w_thresh
                    qcl_mask = qcl[:, height_level].data > qcl_thresh

                    cloud_mask_data[:, h_index, w_index, qcl_index] = (w_mask & qcl_mask).astype(int)

        w_thresh_coord = iris.coords.DimCoord(self.w_threshs, long_name='w_thres', units='m s-1')
        qcl_thresh_coord = iris.coords.DimCoord(self.qcl_threshs, long_name='qcl_thres', units='kg kg-1')
        cloud_mask_cube = iris.cube.Cube(cloud_mask_data, 
                                         long_name='cloud_mask',
                                         dim_coords_and_dims=[(w.coord('time'), 0),
                                                              (w_slice.coord('model_level_number'), 1),
                                                              (w_thresh_coord, 2),
                                                              (qcl_thresh_coord, 3),
                                                              (w.coord('grid_latitude'), 4),
                                                              (w.coord('grid_longitude'), 5)])

        self.results['cloud_mask'] = cloud_mask_cube
