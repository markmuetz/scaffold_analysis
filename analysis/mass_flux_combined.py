import os
from collections import defaultdict

import iris

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MassFluxCombinedAnalysis(Analyzer):
    analysis_name = 'mass_flux_combined'
    multi_file = True

    def set_config(self, config):
	super(MassFluxCombinedAnalysis, self).set_config(config)
	self.model_level_numbers = [int(l) for l in config['model_level_numbers'].split(',')]
	self.w_threshs = [float(t) for t in config['w_threshs'].split(',')]
	self.qcl_threshs = [float(t) for t in config['qcl_threshs'].split(',')]
        self.start_runid = config.getint('start_runid')
	    
    def load(self):
        self.append_log('Override load')
        self.mass_fluxes = defaultdict(list)
        for filename in self.filenames:
            basename = os.path.basename(filename)
            runid = int(basename.split('.')[1])
            if runid >= self.start_runid or True:
                cubes = iris.load(filename)
                for cube in cubes:
                    if cube.name()[:9] == 'mass_flux':
                        (height_level_index, thresh_index) = cube.attributes['mass_flux_key']
                        self.mass_fluxes[(height_level_index, thresh_index)].extend(cube.data)

        self.append_log('Override loaded')

    def run_analysis(self):
        for key, mass_flux in self.mass_fluxes.items():
            (model_level_number, thresh_index) = key

            w_thresh = self.w_threshs[thresh_index]
            qcl_thresh = self.qcl_threshs[thresh_index]
            mf_cube_id = 'mass_flux_z{}_w{}_qcl{}'.format(model_level_number, w_thresh, qcl_thresh)

            values = iris.coords.DimCoord(range(len(mass_flux)), long_name='values')
            mass_flux_cube = iris.cube.Cube(mass_flux, 
                                            long_name=mf_cube_id, 
                                            dim_coords_and_dims=[(values, 0)], 
                                            units='kg s-1')
            mass_flux_cube.attributes['mass_flux_key'] = key
            self.results[mf_cube_id] = mass_flux_cube
