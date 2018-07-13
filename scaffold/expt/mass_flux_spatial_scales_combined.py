import os
from collections import defaultdict
from logging import getLogger

import iris
from omnium import Analyser
from omnium.utils import is_power_of_two

logger = getLogger('scaf.mfssc')


class MassFluxSpatialScalesCombined(Analyser):
    """Combines all mass flux spatial scales

    Only loads files greater than runid == start_runid.
    Loads/saves a separate file for each key (model_level_number, thresh_index, n).
    """
    analysis_name = 'mass_flux_spatial_scales_combined'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.mass_flux_spatial_scales_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.mass_flux_spatial_scales_combined.nc']

    def load(self):
        self.spatial_mass_fluxes = defaultdict(list)
        for filename in self.task.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            if runid >= self.settings.start_runid:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    # N.B. need to deref this to use as a key.
                    (height_index, thresh_index, n) = cube.attributes['mass_flux_spatial_key']
                    key = (height_index, thresh_index, n)
                    assert is_power_of_two(n)

                    self.spatial_mass_fluxes[key].extend(cube.data)
            else:
                logger.debug('skipping runid: {}'.format(runid))

    def run(self):
        for key, spatial_mass_flux in self.spatial_mass_fluxes.items():
            (model_level_number, thresh_index, n) = key

            smf_cube_id = 'spatial-mass-flux_z{0}_w{1}_qcl{1}_n{2}'.format(model_level_number,
                                                                           thresh_index, n)

            values = iris.coords.DimCoord(range(len(spatial_mass_flux)), long_name='values')
            mass_flux_spatial_cube = iris.cube.Cube(spatial_mass_flux,
                                                    long_name=smf_cube_id,
                                                    dim_coords_and_dims=[(values, 0)],
                                                    units='kg s-1')
            mass_flux_spatial_cube.attributes['mass_flux_spatial_key'] = key
            self.results[smf_cube_id] = mass_flux_spatial_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
