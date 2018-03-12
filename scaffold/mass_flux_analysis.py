import matplotlib
import numpy as np

matplotlib.use('Agg')
import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr
from cloud_tracking.utils import label_clds

from scaffold.vertlev import VertLev


class MassFluxAnalyser(Analyser):
    """Works out the mass flux for each cloud.

    Performs some sanity checks on rho/w indices.
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

            # Need to get these to do proper interp of rho onto theta.
            # Work out in 2 ways and check equal because paranoid.
            # AKA sanity checks.
            w_height = w_slice.attributes['heights'][height_level_index]
            w_height2 = vertlevs.z_theta[level_number]

            # N.B. for every 1 w_height, there are 2 rho_heights. Index appropriately.
            rho_height_lower = rho_slice.attributes['heights'][2 * height_level_index]
            rho_height_lower2 = vertlevs.z_rho[level_number - 1]

            rho_height_upper = rho_slice.attributes['heights'][2 * height_level_index + 1]
            rho_height_upper2 = vertlevs.z_rho[level_number]

            # Calc scaling for linear interp.
            alpha = (w_height - rho_height_lower)/(rho_height_upper - rho_height_lower)
            # Paranoia. Well justified it turns out. Saved me from doing wrong analysis.
            assert w_height == w_height2
            assert rho_height_lower == rho_height_lower2
            assert rho_height_upper == rho_height_upper2
            assert 0 <= alpha <= 1

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

                    # Interp rho onto w grid.
                    rho_ss_lower = rho_slice[time_index, height_level_index].data
                    rho_ss_upper = rho_slice[time_index, height_level_index + 1].data

                    rho_ss_interp = (1 - alpha) * rho_ss_lower + alpha * rho_ss_upper

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
