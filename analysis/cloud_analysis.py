import os

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import iris

from omnium.analyzer import Analyzer
from omnium.utils import get_cube
from omnium.consts import Re, L, cp, g


class CloudAnalyzer(Analyzer):
    analysis_name = 'cloud_analysis'

    def set_config(self, config):
	super(CloudAnalyzer, self).set_config(config)
	self.height_levels = [int(l) for l in config['height_levels'].split(',')]
	self.w_threshs = [float(t) for t in config['w_threshs'].split(',')]
	self.qcl_threshs = [float(t) for t in config['qcl_threshs'].split(',')]
	    
    def run_analysis(self):
        cubes = self.cubes

        w = get_cube(cubes, 0, 150)
        qcl = get_cube(cubes, 0, 254)
        rho = get_cube(cubes, 0, 253)

        rho.data = rho.data / Re**2
        rho.units = 'kg m-3'

        rho_heights = []
        for height_level in self.height_levels:
            rho_heights.extend([height_level, height_level + 1])

        w_slice = w[:, self.height_levels].copy()
	self.results['w_slice'] = w_slice
	self.results['qcl_slice'] = qcl[:, self.height_levels].copy()
	self.results['rho_slice'] = rho[:, rho_heights].copy()

        w_thresh_coord = iris.coords.DimCoord(self.w_threshs, long_name='w_thres')
        qcl_thresh_coord = iris.coords.DimCoord(self.qcl_threshs, long_name='qcl_thres')

        cloud_mask_data = np.zeros((w.shape[0], 
                                    len(self.height_levels),
                                    len(self.w_threshs),
                                    len(self.qcl_threshs),
                                    w.shape[-2],
                                    w.shape[-1]))
                              
	for h_index, height_level in enumerate(self.height_levels):
	    for w_index, w_thresh in enumerate(self.w_threshs):
                for qcl_index, qcl_thresh in enumerate(self.qcl_threshs):
                    w_mask = w[:, height_level].data > w_thresh
                    qcl_mask = qcl[:, height_level].data > qcl_thresh

                    cloud_mask_data[:, h_index, w_index, qcl_index] = (w_mask & qcl_mask).astype(int)

        cloud_mask_cube = iris.cube.Cube(cloud_mask_data, 
                                         long_name='cloud_mask',
                                         dim_coords_and_dims=[(w.coord('time'), 0),
                                                              (w_slice.coord('model_level_number'), 1),
                                                              (w_thresh_coord, 2),
                                                              (qcl_thresh_coord, 3),
                                                              (w.coord('grid_latitude'), 4),
                                                              (w.coord('grid_longitude'), 5)])

        self.results['cloud_mask'] = cloud_mask_cube
