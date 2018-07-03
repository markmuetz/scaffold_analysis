import os
from collections import defaultdict
from logging import getLogger

import iris
from omnium.analyser import Analyser
from scaffold.scaffold_settings import settings

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

    settings = settings

    def load(self):
        self.mass_fluxes = defaultdict(list)
        for filename in self.task.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            # if runid >= settings.start_runid:
            if runid >= 24:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:9] == 'mass_flux':
                        (height_level_index, thresh_index) = cube.attributes['mass_flux_key']
                        self.mass_fluxes[(height_level_index, thresh_index)].extend(cube.data)
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
                                            units='kg s-1')
            mass_flux_cube.attributes['mass_flux_key'] = key
            self.results[mf_cube_id] = mass_flux_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
