import os
from collections import defaultdict
from logging import getLogger

import iris
from omnium.analyzer import Analyzer

logger = getLogger('om.oc')

class OrgCombined(Analyzer):
    analysis_name = 'org_combined'
    multi_file = True

    def set_config(self, config):
        super(OrgCombined, self).set_config(config)
        self.start_runid = config.getint('start_runid')
            
    def load(self):
        self.append_log('Override load')
        self.dists = defaultdict(list)
        for filename in self.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            if runid >= self.start_runid:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:4] == 'dist':
                        (height_level_index, thresh_index) = cube.attributes['dist_key']
                        self.dists[(height_level_index, thresh_index)].extend(cube.data)
            else:
                logger.debug('skipping runid: {}'.format(runid))

        self.append_log('Override loaded')

    def run_analysis(self):
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
