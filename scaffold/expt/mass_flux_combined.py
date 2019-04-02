import os
from collections import defaultdict
from logging import getLogger

import iris

from omnium import Analyser

logger = getLogger('scaf.mfc')

class MassFluxCombinedAnalysis(Analyser):
    """Combines previously worked out mass_fluxes into one large sequence of mass_fluxes.

    load(...) only loads mass fluxes greater than a given start_runid.
    """
    analysis_name = 'mass_flux_combined'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.mass_flux_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.mass_flux_combined.nc']

    def load(self):
        self.mass_fluxes = defaultdict(list)
        self.total_mass_fluxes = defaultdict(list)
        self.sigmas = defaultdict(list)
        for filename in self.task.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            # if self.settings.start_runid <= runid < self.settings.end_runid:
            if 0 <= runid < 76:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:9] == 'mass_flux':
                        (height_level_index, thresh_index) = cube.attributes['mass_flux_key']
                        self.mass_fluxes[(height_level_index, thresh_index)].extend(cube.data)
                    elif cube.name().startswith('total_mass_flux'):
                        (height_level_index, thresh_index) = cube.attributes['total_mass_flux_key']
                        self.total_mass_fluxes[(height_level_index, thresh_index)].extend(cube.data)
                    elif cube.name().startswith('sigma'):
                        (height_level_index, thresh_index) = cube.attributes['sigma_key']
                        self.sigmas[(height_level_index, thresh_index)].extend(cube.data)
            else:
                logger.debug('skipping runid: {}'.format(runid))

    def run(self):
        for key, mass_flux in self.mass_fluxes.items():
            (model_level_number, thresh_index) = key

            mf_cube_id = 'mass_flux_z{0}_w{1}_qcl{1}'.format(model_level_number, thresh_index)

            values = iris.coords.DimCoord(range(len(mass_flux)), long_name='values')
            mass_flux_cube = iris.cube.Cube(mass_flux,
                                            long_name=mf_cube_id,
                                            dim_coords_and_dims=[(values, 0)],
                                            units='kg m-2 s-1')
            mass_flux_cube.attributes['mass_flux_key'] = key
            self.results[mf_cube_id] = mass_flux_cube

        for key, total_mass_flux in self.total_mass_fluxes.items():
            (model_level_number, thresh_index) = key

            total_mf_cube_id = 'total_mass_flux_z{0}_w{1}_qcl{1}'.format(model_level_number, thresh_index)

            values = iris.coords.DimCoord(range(len(total_mass_flux)), long_name='values')
            total_mass_flux_cube = iris.cube.Cube(total_mass_flux,
                                                  long_name=total_mf_cube_id,
                                                  dim_coords_and_dims=[(values, 0)],
                                                  units='kg m-2 s-1')
            total_mass_flux_cube.attributes['total_mass_flux_key'] = key
            self.results[total_mf_cube_id] = total_mass_flux_cube

        for key, sigma in self.sigmas.items():
            (model_level_number, thresh_index) = key

            sigma_cube_id = 'sigma_z{0}_w{1}_qcl{1}'.format(model_level_number, thresh_index)

            values = iris.coords.DimCoord(range(len(sigma)), long_name='values')
            sigma_cube = iris.cube.Cube(sigma,
                                        long_name=sigma_cube_id,
                                        dim_coords_and_dims=[(values, 0)],
                                        units='')
            sigma_cube.attributes['sigma_key'] = key
            self.results[sigma_cube_id] = sigma_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
