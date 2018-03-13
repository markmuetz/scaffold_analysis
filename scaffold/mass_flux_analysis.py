import matplotlib
import numpy as np

matplotlib.use('Agg')
import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr

from scaffold.vertlev import VertLev
from scaffold.utils import interp_vert_rho2w


class MassFluxAnalyser(Analyser):
    """Works out the mass flux for each cloud.

    Uses the output from cloud_analysis, and labels each of the contiguous clouds
    (with diagonal=True). Uses these clouds to work out the mass flux for each cloud.

    """
    analysis_name = 'mass_flux_analysis'
    single_file = True

    def run_analysis(self):
        cubes = self.cubes

        w_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'w_slice')
        rho_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'rho_slice')

        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        qcl_thresh_coord = cloud_mask_cube.coord('qcl_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')

        vertlevs = VertLev(self.suite.suite_dir)
        # level_number refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, level_number in enumerate(level_number_coord.points):
            for thresh_index in range(w_thresh_coord.shape[0]):
                # N.B. I just take the diagonal indices.
                w_thresh = w_thresh_coord.points[thresh_index]
                qcl_thresh = qcl_thresh_coord.points[thresh_index]
                labelled_clouds_cube_id = 'labelled_clouds_z{}_w{}_qcl{}'.format(level_number,
                                                                                 w_thresh,
                                                                                 qcl_thresh)
                labelled_clouds_cube = get_cube_from_attr(cubes, 'omnium_cube_id', labelled_clouds_cube_id)

                mass_fluxes = []
                for time_index in range(cloud_mask_cube.data.shape[0]):
                    # w_ss == w_snapshot.
                    w_ss = w_slice[time_index, height_level_index].data

                    rho_ss_interp = interp_vert_rho2w(vertlevs, w_slice, rho_slice, time_index,
                                                      height_level_index, level_number)

                    labelled_clouds_ss = labelled_clouds_cube[time_index].data.astype(int)
                    mf_ss = rho_ss_interp * w_ss
                    max_labelled_cloud_index = np.max(labelled_clouds_ss)

                    for i in range(1, max_labelled_cloud_index + 1):
                        mask = (labelled_clouds_ss == i)
                        mass_flux = mf_ss[mask].sum()
                        mass_fluxes.append(mass_flux)

                mf_cube_id = 'mass_flux_z{}_w{}_qcl{}'.format(level_number, w_thresh, qcl_thresh)

                values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
                mass_flux_cube = iris.cube.Cube(mass_fluxes,
                                                long_name=mf_cube_id,
                                                dim_coords_and_dims=[(values, 0)],
                                                units='kg s-1')

                mass_flux_cube.attributes['mass_flux_key'] = (height_level_index, thresh_index)
                self.results[mf_cube_id] = mass_flux_cube

