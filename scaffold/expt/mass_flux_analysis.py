import matplotlib
import numpy as np

matplotlib.use('Agg')
import iris

from omnium import Analyser
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
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.cloud_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.{runid:03}.mass_flux_analysis.nc']
    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).cloud_analysis.nc'

    def load(self):
        self.load_cubes()

    def run(self):
        cubes = self.cubes

        w_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'w_slice')
        rho_slice = get_cube_from_attr(cubes, 'omnium_cube_id', 'rho_slice')

        cloud_mask_id = 'cloud_mask'
        cloud_mask_cube = get_cube_from_attr(cubes, 'omnium_cube_id', cloud_mask_id)

        w_thresh_coord = cloud_mask_cube.coord('w_thres')
        qcl_thresh_coord = cloud_mask_cube.coord('qcl_thres')
        level_number_coord = cloud_mask_cube.coord('model_level_number')

        num_domain_grid_cells = w_slice.shape[2] * w_slice.shape[3]

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
                total_mass_fluxes = []
                sigmas = []

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
                    total_mass_fluxes.append(mf_ss[labelled_clouds_ss >= 1].sum())
                    sigmas.append((labelled_clouds_ss >= 1).sum() / num_domain_grid_cells)

                mf_cube_id = 'mass_flux_z{}_w{}_qcl{}'.format(level_number, w_thresh, qcl_thresh)
                values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
                mass_flux_cube = iris.cube.Cube(mass_fluxes,
                                                long_name=mf_cube_id,
                                                dim_coords_and_dims=[(values, 0)],
                                                units='kg m-2 s-1')

                total_mf_cube_id = 'total_mass_flux_z{}_w{}_qcl{}'.format(level_number, w_thresh, qcl_thresh)
                values = iris.coords.DimCoord(range(len(total_mass_fluxes)), long_name='values')
                total_mass_flux_cube = iris.cube.Cube(total_mass_fluxes,
                                                      long_name=total_mf_cube_id,
                                                      dim_coords_and_dims=[(values, 0)],
                                                      units='kg m-2 s-1')

                mass_flux_cube.attributes['mass_flux_key'] = (height_level_index, thresh_index)
                total_mass_flux_cube.attributes['total_mass_flux_key'] = (height_level_index, thresh_index)

                sigma_cube_id = 'sigma_z{}_w{}_qcl{}'.format(level_number, w_thresh, qcl_thresh)
                values = iris.coords.DimCoord(range(len(sigmas)), long_name='values')
                sigma_cube = iris.cube.Cube(sigmas,
                                            long_name=sigma_cube_id,
                                            dim_coords_and_dims=[(values, 0)],
                                            units='')

                mass_flux_cube.attributes['mass_flux_key'] = (height_level_index, thresh_index)
                total_mass_flux_cube.attributes['total_mass_flux_key'] = (height_level_index, thresh_index)
                sigma_cube.attributes['sigma_key'] = (height_level_index, thresh_index)

                self.results[mf_cube_id] = mass_flux_cube
                self.results[total_mf_cube_id] = total_mass_flux_cube
                self.results[sigma_cube_id] = sigma_cube

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
