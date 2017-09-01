import os
from collections import defaultdict
from logging import getLogger

import iris
from omnium.analyser import Analyser

logger = getLogger('om.mfc')

class MassFluxCombinedAnalysis(Analyser):
    analysis_name = 'mass_flux_combined'
    multi_file = True

    def set_config(self, config):
        super(MassFluxCombinedAnalysis, self).set_config(config)
        self.start_runid = config.getint('start_runid')
            
    def load(self):
        self.append_log('Override load')
        self.mass_fluxes = defaultdict(list)
        for filename in self.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            if runid >= self.start_runid:
                logger.debug('adding runid: {}'.format(runid))
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:9] == 'mass_flux':
                        (height_level_index, thresh_index) = cube.attributes['mass_flux_key']
                        self.mass_fluxes[(height_level_index, thresh_index)].extend(cube.data)
            else:
                logger.debug('skipping runid: {}'.format(runid))

        self.append_log('Override loaded')

    def run_analysis(self):
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

