import itertools
from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser
from omnium.utils import get_cube
from omnium.consts import L, cp, Lf

logger = getLogger('scaf.dump_entr')

SCALARS = ['MSE', 'FMSE']
CLD_OPTS = ['scaffold', 'scaffold_plus', 'scaffold_plus_strong', 'SS2012', 'swann_ud', 'swann_core']
METHODS = ['normal', 'acc']

class DumpEntr(Analyser):
    analysis_name = 'dump_entr'
    multi_expt = True

    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da480.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.dump_entr.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self._plot_entrainment()
        plt.close('all')

    def _plot_entrainment(self):
        opts = []
        for p in itertools.product(SCALARS, CLD_OPTS, METHODS):
            opt = {'SCALAR': p[0],
                    'CLD_DEF': p[1],
                    'METHOD': p[2]}
            opts.append(opt)

        for i, expt in enumerate(self.task.expts):
            da = self.expt_cubes[expt]

            theta_cube = get_cube(da, 0, 4)
            q_cube = get_cube(da, 0, 10)

            qcl_cube = get_cube(da, 0, 254)
            qcf_cube = get_cube(da, 0, 12)
            qcf2_cube = get_cube(da, 0, 271)
            qgr_cube = get_cube(da, 0, 273)
            w_cube = get_cube(da, 0, 150)

            z = theta_cube.coord('atmosphere_hybrid_height_coordinate').points

            # Calc MSE (h), FMSE (h2)
            theta = theta_cube.data
            q = q_cube.data
            qcf = qcf_cube.data
            qcf2 = qcf2_cube.data
            qgr = qgr_cube.data
            h = theta + L / cp * q
            h2 = theta + L / cp * q - Lf / cp * (qcf + qcf2 + qgr)

            w = w_cube.data
            qcl = qcl_cube.data
            qcf = qcf_cube.data
            qcf2 = qcf2_cube.data
            theta = theta_cube.data

            for opt in opts:
                logger.debug('entr for {}: {}', expt, opt)
                if opt['CLD_DEF'] == 'SS2012':
                    # Using def from Stirling and Stratton 2012:
                    # qcl or qcf (or qcf2) > 1e-5 kg/kg
                    # w > 0
                    # +ve buoyancy (here defined as theta > theta.mean())
                    qcl_thresh = 1e-5
                    w_thresh = 0
                    env_mask = (((qcl > qcl_thresh) | (qcf > qcl_thresh) | (qcf2 > qcl_thresh))
                                & (w > w_thresh)
                                & (theta > theta.mean(axis=(1, 2))[:, None, None]))
                elif opt['CLD_DEF'] == 'swann_ud':
                    qcl_thresh = 1e-6
                    w_thresh = 0
                    env_mask = (((qcl > qcl_thresh) | (qcf > qcl_thresh) | (qcf2 > qcl_thresh))
                                & (w > w_thresh))
                elif opt['CLD_DEF'] == 'swann_core':
                    qcl_thresh = 1e-6
                    w_thresh = 0
                    env_mask = (((qcl > qcl_thresh) | (qcf > qcl_thresh) | (qcf2 > qcl_thresh))
                                & (w > w_thresh)
                                & (theta > theta.mean(axis=(1, 2))[:, None, None]))
                elif opt['CLD_DEF'] == 'scaffold':
                    qcl_thresh = 5e-6
                    w_thresh = 1
                    env_mask = (qcl > qcl_thresh) & (w > w_thresh)
                elif opt['CLD_DEF'] == 'scaffold_plus':
                    qcl_thresh = 5e-6
                    w_thresh = 1
                    env_mask = (((qcl > qcl_thresh) | (qcf > qcl_thresh) | (qcf2 > qcl_thresh))
                                & (w > w_thresh))
                elif opt['CLD_DEF'] == 'scaffold_plus_strong':
                    qcl_thresh = 5e-6
                    w_thresh = 5
                    env_mask = (((qcl > qcl_thresh) | (qcf > qcl_thresh) | (qcf2 > qcl_thresh))
                                & (w > w_thresh))
                cld_mask = ~env_mask

                # Calc entr.
                if opt['SCALAR'] == 'MSE':
                    scalar = h
                elif opt['SCALAR'] == 'FMSE':
                    scalar = h2
                # Can't do without knowing source terms of qt.
                # elif opt['SCALAR'] == 'qt':
                #     scalar = (expt_res['q_cube'].data +
                #         expt_res['qcl_cube'].data +
                #         expt_res['qcf_cube'].data +
                #         expt_res['qcf2_cube'].data)

                if opt['METHOD'] == 'acc':
                    scalar_c = (np.ma.masked_array(scalar * w, cld_mask).mean(axis=(1, 2))
                                / np.ma.masked_array(w, cld_mask).mean(axis=(1, 2)))
                    scalar_e = np.ma.masked_array(scalar, env_mask).mean(axis=(1, 2))

                    dscalar_c_dz = (scalar_c[2:] - scalar_c[:-2])/(z[2:] - z[:-2])
                    entr = dscalar_c_dz / (scalar_e[1:-1] - scalar_c[1:-1])

                elif opt['METHOD'] == 'normal':
                    scalar_c = np.ma.masked_array(scalar, cld_mask).mean(axis=(1, 2))
                    scalar_e = np.ma.masked_array(scalar, env_mask).mean(axis=(1, 2))

                    dscalar_c_dz = (scalar_c[2:] - scalar_c[:-2]) / (z[2:] - z[:-2])
                    entr = dscalar_c_dz / (scalar_e[1:-1] - scalar_c[1:-1])

                # Plot results
                z_km = z / 1000

                figname = 'entr_{}_{}_{}_{}'.format(expt,
                                                    opt['SCALAR'],
                                                    opt['CLD_DEF'],
                                                    opt['METHOD'])
                fig = plt.figure(figname)
                plt.clf()
                fig, axgrid = plt.subplots(1, 5, fig=fig, sharey=True, num=figname, figsize=(15, 12))
                ax0, ax1, ax2, ax3, ax4 = axgrid

                # ~10 km
                # ztop_index = 55
                ztop_index = 75
                # full.
                # ztop_index = 98
                ax0.plot(scalar_c[1:ztop_index], z_km[1:ztop_index], label='cld')
                ax0.plot(scalar_e[1:ztop_index], z_km[1:ztop_index], label='env')
                ax0.legend()

                ax1.plot(dscalar_c_dz[:ztop_index - 1] * 1000, z_km[1:ztop_index])

                ax2.plot(scalar_e[1:ztop_index] - scalar_c[1:ztop_index], z_km[1:ztop_index])

                ax3.plot(env_mask.sum(axis=(1, 2)) / (env_mask.shape[1] * env_mask.shape[2]), z_km)
                ax3.set_xlim((0, 0.1))

                ax4.plot(entr[:ztop_index - 1] * 1000, z_km[1:ztop_index])
                ax4.set_xlim((-1, 1))
                ax4.axvline(x=0, color='k', linestyle='--')

                ax0.set_ylabel('height (km)')
                ax0.set_ylim((0, 15))
                ax0.set_xlabel('{}'.format(opt['SCALAR']))
                ax1.set_xlabel('$\\frac{{d{}_c}}{{dz}}$'.format(opt['SCALAR']))
                ax2.set_xlabel('${0}_e - {0}_c$'.format(opt['SCALAR']))
                ax3.set_xlabel('$\\sigma$')
                ax4.set_xlabel('$\\epsilon$ (km$^{-1}$)')

                plt.savefig(self.file_path('entr_' + figname + '.png'))
