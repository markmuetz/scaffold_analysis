import os
from collections import defaultdict
from logging import getLogger

import iris

from omnium import Analyser

logger = getLogger('scaf.oc')


class OrgCombined(Analyser):
    """Combines the organizational data for runid >= self.start_runid."""
    analysis_name = 'org_combined'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.org_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.org_combined.nc']
    # TODO: runid

    def load(self):
        self.dists = defaultdict(list)
        for filename in self.task.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            if runid >= self.settings.start_runid:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:4] == 'dist':
                        (height_level_index, thresh_index) = cube.attributes['dist_key']
                        self.dists[(height_level_index, thresh_index)].extend(cube.data)
            else:
                logger.debug('skipping runid: {}'.format(runid))

    def run(self):
        for key, dists in self.dists.items():
            (model_level_number, thresh_index) = key

            dist_cube_id = 'dist_{0}_w{1}_qcl{1}'.format(model_level_number, thresh_index)

            values = iris.coords.DimCoord(range(len(dists)), long_name='values')
            dist_cube = iris.cube.Cube(dists,
                                       long_name=dist_cube_id,
                                       dim_coords_and_dims=[(values, 0)],
                                       units='m')
            dist_cube.attributes['dist_key'] = key
            self.results[dist_cube_id] = dist_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
