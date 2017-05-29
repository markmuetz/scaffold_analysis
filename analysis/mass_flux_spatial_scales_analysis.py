import os
from collections import OrderedDict
from logging import getLogger

import numpy as np
import iris

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr, coarse_grain

logger = getLogger('omnium')


class MassFluxSpatialScalesAnalyzer(Analyzer):
    analysis_name = 'mass_flux_spatial_scales_analysis'

    def run_analysis(self):
        cubes = self.cubes

        w = get_cube_from_attr(cubes, 'omnium_cube_id', 'w_slice')
	cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', 'cloud_mask')

	mass_flux = OrderedDict()
	for time_index in range(cloud_mask_cube.shape[0]):
            logger.debug('Time: {}/{}'.format(time_index + 1, cloud_mask_cube.shape[0]))

            for height_index in range(cloud_mask_cube.shape[1]):
                # One value of w (w SnapShot) for each time/height.
                w_ss = w[time_index, height_index].data

                # There are 3 values of cloud_mask for each time/height:
                for thresh_index in range(cloud_mask_cube.shape[2]):
                    # N.B. take diagonal of thresh, i.e. low/low, med/med, hi/hi.
                    cloud_mask_ss = cloud_mask_cube[time_index, height_index, 
                                                    thresh_index, thresh_index].data.astype(bool)
                    coarse_data = coarse_grain(w_ss, cloud_mask_ss)
                    for n, coarse_datum in coarse_data:
                        key = (height_index, thresh_index, n)
                        # logger.debug('Storing coarse data for: {}'.format(key))
                        if key not in mass_flux:
                            mass_flux[key] = []
                        else:
                            mass_flux[key] = np.append(mass_flux[key], coarse_datum.flatten())
	
	for key, mass_fluxes in mass_flux.items():
            logger.debug('Building iris cube for: {}'.format(key))
	    values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
	    name = 'spatial-mass-flux-h{0}_n{1}_w{2}_qcl{2}'.format(*key)
	    mass_flux_cube = iris.cube.Cube(mass_fluxes, 
					    long_name=name, 
					    dim_coords_and_dims=[(values, 0)], 
					    units='kg s-1')
	    self.results[name] = mass_flux_cube
