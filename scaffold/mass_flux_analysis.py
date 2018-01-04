import matplotlib
import numpy as np

matplotlib.use('Agg')
import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr
from cloud_tracking.utils import label_clds

from scaffold.vertlev import VertLev


class MassFluxAnalyser(Analyser):
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
        # height_level refers to orig cube.
        # height_level_index refers to w as it has already picked out the height levels.
        for height_level_index, height_level in enumerate(level_number_coord.points):
            # Need to get these to do proper interp of rho onto theta.
            # Work out in 2 ways and check equal because paranoid.
            w_height = w_slice.attributes['heights'][height_level_index]
            w_height2 = vertlevs.z_theta[height_level]

            # N.B. for every 1 w_height, there are 2 rho_heights. Index appropriately.
            rho_height_lower = rho_slice.attributes['heights'][2 * height_level_index]
            rho_height_lower2 = vertlevs.z_rho[height_level - 1]

            rho_height_upper = rho_slice.attributes['heights'][2 * height_level_index + 1]
            rho_height_upper2 = vertlevs.z_rho[height_level]

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

                blob_cube = w_slice[:, height_level_index].copy()
                blob_cube.units = ''

                blob_cube_data = np.zeros_like(blob_cube.data)
                mass_fluxes = []
                for time_index in range(cloud_mask_cube.data.shape[0]):
                    w_ss = w_slice[time_index, height_level_index].data
                    rho_ss_lower = rho_slice[time_index, height_level_index].data
                    rho_ss_upper = rho_slice[time_index, height_level_index + 1].data

                    # rho_ss_interp_old = (rho_ss_lower + rho_ss_upper) / 2
                    rho_ss_interp = (1 - alpha) * rho_ss_lower + alpha* rho_ss_upper
                    #rho_ss_interp = (rho_ss_lower + rho_ss_upper) / 2
                    cloud_mask_ss = cloud_mask_cube[time_index,
                                                    height_level_index,
                                                    thresh_index,
                                                    thresh_index].data.astype(bool)
                    max_blob_index, blobs = label_clds(cloud_mask_ss, True)
                    blob_cube_data[time_index] = blobs
                    mf_ss = rho_ss_interp * w_ss

                    for i in range(1, max_blob_index + 1):
                        mask = (blobs == i)
                        mass_flux = mf_ss[mask].sum()
                        mass_fluxes.append(mass_flux)

                mf_cube_id = 'mass_flux_z{}_w{}_qcl{}'.format(height_level, w_thresh, qcl_thresh)
                blob_cube_id = 'blob_z{}_w{}_qcl{}'.format(height_level, w_thresh, qcl_thresh)

                blob_cube.rename(blob_cube_id)
                blob_cube.data = blob_cube_data
                self.results[blob_cube_id] = blob_cube

                values = iris.coords.DimCoord(range(len(mass_fluxes)), long_name='values')
                mass_flux_cube = iris.cube.Cube(mass_fluxes, 
                                                long_name=mf_cube_id, 
                                                dim_coords_and_dims=[(values, 0)], 
                                                units='kg s-1')

                mass_flux_cube.attributes['mass_flux_key'] = (height_level_index, thresh_index)
                self.results[mf_cube_id] = mass_flux_cube
