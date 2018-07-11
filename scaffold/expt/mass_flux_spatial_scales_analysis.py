from collections import OrderedDict
from logging import getLogger

import iris
import numpy as np

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr, coarse_grain

from scaffold.vertlev import VertLev
from scaffold.utils import interp_vert_rho2w

logger = getLogger('scaf.mfssa')


class MassFluxSpatialScalesAnalyser(Analyser):
    """Split mass flux into powers of 2 subdomains, calc the total conv. mf for each subdomain.

    conv. mf = mass flux where cloud_mask == true.
    Saves a separate file for each key (model_level_number, thresh_index, n).
    """
    analysis_name = 'mass_flux_spatial_scales_analysis'
    single_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.cloud_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.{runid:03}.mass_flux_spatial_scales_analysis.nc']

    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).cloud_analysis.nc'

    def load(self):
        self.load_cubes()

    def run(self):
        cubes = self.cubes

        w_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'w_slice')
        rho_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'rho_slice')

        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', 'cloud_mask')

        level_number_coord = cloud_mask_cube.coord('model_level_number')
        vertlevs = VertLev(self.suite.suite_dir)

        mass_flux = OrderedDict()
        for time_index in range(cloud_mask_cube.shape[0]):
            logger.debug('time: {}/{}'.format(time_index + 1, cloud_mask_cube.shape[0]))

            for height_level_index, level_number in enumerate(level_number_coord.points):
                # One value of w_slice (w SnapShot) for each time/height.
                w_ss = w_slice[time_index, height_level_index].data

                rho_ss_interp = interp_vert_rho2w(vertlevs, w_slice, rho_slice, time_index,
                                                  height_level_index, level_number)

                mf_ss = rho_ss_interp * w_ss
                N = mf_ss.shape[-1]

                # There are 3 values of cloud_mask for each time/height:
                for thresh_index in range(cloud_mask_cube.shape[2]):
                    # N.B. take diagonal of thresh, i.e. low/low, med/med, hi/hi.
                    cloud_mask_ss = cloud_mask_cube[time_index, height_level_index,
                                                    thresh_index, thresh_index].data.astype(bool)
                    # Heart of the analysis. Coarse grain the data.
                    # i.e. split into 4, then 8, then 16... subdomains.
                    # For each subdomain (and each power), save the convective mass flux for that
                    # subdomain (convective == where cloud_mask is true). Store all of these.
                    coarse_data = coarse_grain(mf_ss, cloud_mask_ss, self.settings.npow)
                    for n, coarse_datum in coarse_data:
                        key = (height_level_index, thresh_index, n)
                        Nsubdom = N / n
                        if key not in mass_flux:
                            mass_flux[key] = []
                        mass_flux[key] = np.append(mass_flux[key], coarse_datum.flatten() / Nsubdom**2)

        for key, mass_fluxes in mass_flux.items():
            logger.debug('building iris cube for: {}'.format(key))
            values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
            name = 'spatial-mass-flux-h{0}_w{1}_qcl{1}_n{2}'.format(*key)
            mass_flux_spatial_cube = iris.cube.Cube(mass_fluxes,
                                                    long_name=name,
                                                    dim_coords_and_dims=[(values, 0)],
                                                    units='kg s-1')

            mass_flux_spatial_cube.attributes['mass_flux_spatial_key'] = key
            self.results[name] = mass_flux_spatial_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
