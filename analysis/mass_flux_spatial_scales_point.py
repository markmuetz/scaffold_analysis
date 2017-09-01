from collections import OrderedDict
from logging import getLogger

import iris
import numpy as np
from cloud_tracking.utils import label_clds
from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr

logger = getLogger('om.mfssp')


class MassFluxSpatialScalesPoint(Analyser):
    analysis_name = 'mass_flux_spatial_scales_point'

    def run_analysis(self):
        cubes = self.cubes

        w = get_cube_from_attr(cubes, 'omnium_cube_id', 'w_slice')
        rho_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'rho_slice')
        lat = w.coord('grid_latitude').points
        lon = w.coord('grid_longitude').points

        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', 'cloud_mask')

        mass_flux = OrderedDict()
        for time_index in range(cloud_mask_cube.shape[0]):
            logger.debug('Time: {}/{}'.format(time_index + 1, cloud_mask_cube.shape[0]))

            for height_index in range(cloud_mask_cube.shape[1]):
                # One value of w (w SnapShot) for each time/height.
                w_ss = w[time_index, height_index].data

                rho_ss_lower = rho_slice[time_index, height_index].data
                rho_ss_upper = rho_slice[time_index, height_index + 1].data
                rho_ss_interp = (rho_ss_lower + rho_ss_upper) / 2

                mf_ss = rho_ss_interp * w_ss

                # There are 3 values of cloud_mask for each time/height:
                for thresh_index in range(cloud_mask_cube.shape[2]):
                    # N.B. take diagonal of thresh, i.e. low/low, med/med, hi/hi.
                    cloud_mask_ss = cloud_mask_cube[time_index, height_index, 
                                                    thresh_index, thresh_index].data.astype(bool)

                    max_blob_index, blobs = label_clds(cloud_mask_ss, True)
                    cloud_coords = []
                    for i in range(1, max_blob_index + 1):
                        blob_indices = np.where(b == i)
                        cloud_lat = lat[blob_indices[0]].mean()
                        cloud_lon = lon[blob_indices[1]].mean()
                        cloud_coords.append((cloud_lat, cloud_lon))

                    lat_values = iris.coords.DimCoord(range(cloud_coords.shape[0]), 
                                                      long_name='lat_values')
                    lon_values = iris.coords.DimCoord(range(cloud_coords.shape[0]), 
                                                      long_name='lon_values')

                    lat_cube_id = ''
                    lat_cube = iris.cube.Cube(cloud_coords[:, 0], 
                                              long_name=lat_cube_id, 
                                              dim_coords_and_dims=[(lat_values, 0)], 
                                              units='kg s-1')

                    

        
        for key, mass_fluxes in mass_flux.items():
            logger.debug('Building iris cube for: {}'.format(key))
            values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
            name = 'spatial-mass-flux-h{0}_w{1}_qcl{1}_n{2}'.format(*key)
            mass_flux_spatial_cube = iris.cube.Cube(mass_fluxes, 
                                            long_name=name, 
                                            dim_coords_and_dims=[(values, 0)], 
                                            units='kg s-1')

            mass_flux_spatial_cube.attributes['mass_flux_spatial_key'] = key
            self.results[name] = mass_flux_spatial_cube
