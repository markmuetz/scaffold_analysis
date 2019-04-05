import os
from glob import glob
from logging import getLogger

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import iris

from omnium import Analyser
from omnium.utils import cd, get_cube_from_attr, cm_to_inch


logger = getLogger('scaf.summ_corr')

SORTED_EXPTS = [
    'S0W0Forced',
    'S4W0Forced',
    'S0W5Forced',
    'S4W5Forced',
    'RWP_C1',
    'RWP_C2',
    'RWP_C3',
    'RWP_C4',
    'RWP_C5',
    'RWP_C6',
    'RWP_C7',
    'RWP_C8',
    'RWP_C9',
    'RWP_C10',
]


def get_mass_flux(expt):
    mfc = iris.load('{}/atmos.mass_flux_combined.nc'.format(expt))
    mf = get_cube_from_attr(mfc, 'omnium_cube_id', 'mass_flux_z1_w1_qcl1')
    tmf = get_cube_from_attr(mfc, 'omnium_cube_id', 'total_mass_flux_z1_w1_qcl1')
    sigma = get_cube_from_attr(mfc, 'omnium_cube_id', 'sigma_z1_w1_qcl1')
    return {'num_clds': mf.shape[0], 'mfpc': mf.data.mean(),
            'tmf': tmf.data.mean(), 'sigma': sigma.data.mean()}


def get_cloud_lifetimes(expt):
    with open('{}/atmos.cloud_tracking_z1_t1.txt'.format(expt), 'r') as f:
        lines = f.readlines()
    # print(expt)
    # print(lines)

    simple_lifetime = float(lines[4].split(',')[-1].strip())
    complex_lifetime = float(lines[7].split(',')[-1].strip())
    return {'simple_lifetime': simple_lifetime, 'complex_lifetime': complex_lifetime}


def savefig(cwd):
    plt.savefig('{}/summary_info_{}.pdf'.format(cwd, plt.gcf().get_label()))


class SummaryCorrelations(Analyser):
    analysis_name = 'summary_correlations'
    single_file = True

    # Completely ignore omnium file loading mechanism!
    input_dir = 'omnium_output/{version_dir}'
    input_filename = '{input_dir}/summary_correlations.dummy'
    output_dir = 'omnium_output/{version_dir}/suite_summary_correlations'
    output_filenames = ['{output_dir}/atmos.summary_correlations.dummy']

    def load(self):
        self.basedir = os.path.dirname(self.task.filenames[0])
        self.outputdir = os.path.dirname(self.task.output_filenames[0])

        with cd(self.basedir):
            si_fns = sorted(glob('*/atmos.restart_dump_summary_info.hdf'))
            self.df_org = pd.read_csv('suite_S0W0Forced_S4W0Forced_S0W5Forced_S4W5Forced_'
                                 'RWP_C1_RWP_C2_RWP_C3_RWP_C4_RWP_C5_'
                                 'RWP_C6_RWP_C7_RWP_C8_RWP_C9_RWP_C10'
                                 '/org_data.csv')
            self.df_tds = {}
            self.df_dyns = {}
            for fn in si_fns:
                expt = os.path.dirname(fn)
                self.df_tds[expt] = pd.read_hdf(fn, 'thermodynamic_summary')
                self.df_dyns[expt] = pd.read_hdf(fn, 'dynamic_summary')

            self.cloud_lifetimes = {}
            self.mass_flux = {}
            for expt in SORTED_EXPTS:
                self.cloud_lifetimes[expt] = get_cloud_lifetimes(expt)
                self.mass_flux[expt] = get_mass_flux(expt)

    def run(self):
        for expt in SORTED_EXPTS:
            print(expt)
            df_td = self.df_tds[expt]
            # Skip first column: runid.
            for col in df_td.columns[1:]:
                print('  {}: {}, {}'.format(col, df_td[col].mean(), df_td[col].std()))
            df_dyn = self.df_dyns[expt]
            # Skip first column: runid.
            for col in df_dyn.columns[1:]:
                print('  {}: {}, {}'.format(col, df_dyn[col].mean(), df_dyn[col].std()))

        cols = (list(df_dyn.columns[1:]) + list(df_td.columns[1:]) +
                list(self.df_org.columns[2:]) +
                ['num_clds', 'mfpc', 'tmf', 'sigma'] +
                ['simple_lifetime', 'complex_lifetime'])

        all_col_values = {}
        for i, col in enumerate(cols):
            # plt.figure(col)
            # plt.clf()
            col_values = []
            for expt in SORTED_EXPTS:
                if col in df_dyn.columns[1:]:
                    col_values.append(self.df_dyns[expt][col].mean())
                if col in df_td.columns[1:]:
                    col_values.append(self.df_tds[expt][col].mean())
                elif col in self.df_org.columns[2:]:
                    col_values.append(self.df_org[(self.df_org.expt == expt) &
                                             (self.df_org.group == 1)][col].values[0])
                elif col in ['num_clds', 'mfpc', 'tmf', 'sigma']:
                    col_values.append(self.mass_flux[expt][col])
                elif col in ['simple_lifetime', 'complex_lifetime']:
                    col_values.append(self.cloud_lifetimes[expt][col])

            # plt.plot(SORTED_EXPTS, col_values, 'ko')
            # if col in ['LCL', 'LFC', 'LNB']:
            #     plt.gca().invert_yaxis()
            # savefig(cwd)
            all_col_values[col] = col_values

        corr = np.zeros((len(cols), len(cols)))
        pvals = np.zeros((len(cols), len(cols)))
        # f, axes = plt.subplots(len(cols), len(cols), num='corr_plots')

        for i, col1 in enumerate(cols):
            for j, col2 in enumerate(cols):
                d1 = np.array(all_col_values[col1])
                d2 = np.array(all_col_values[col2])
                if True:
                    if col1 in ['LCL', 'LFC', 'LNB']:
                        d1 = -d1
                    if col2 in ['LCL', 'LFC', 'LNB']:
                        d2 = -d2
                    if col1 in ['CIN']:
                        d1 = np.abs(d1)
                    if col2 in ['CIN']:
                        d2 = np.abs(d2)

                lr = stats.linregress(d1, d2)
                corr[i, j] = lr.rvalue
                pvals[i, j] = lr.pvalue
                # axes[i, j].plot(d1, d2, 'kx')
        # savefig(cwd)

        self.cols = cols
        self.all_col_values = all_col_values
        self.corr = corr
        self.pvals = pvals

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def _plot_matrix(self, name, cols, data, cmap=None, sum_rows=False):
        plt.figure(name, figsize=cm_to_inch(29.7, 21.0))
        plt.clf()
        if cmap:
            plt.imshow(data, cmap=cmap)
        else:
            plt.imshow(data)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                plt.annotate('{:.2f}'.format(data[i, j]), (i - 0.3, j + 0.1),
                             color='w', fontsize=8)

        if sum_rows:
            for j in range(data.shape[1]):
                plt.annotate('{:.2f}'.format(data[:, j].sum()), (data.shape[0] - 0.3, j + 0.1),
                             color='k', fontsize=8, annotation_clip=False)

        plt.xticks(range(len(cols)), cols, rotation=90)
        plt.yticks(range(len(cols)), cols)
        plt.tight_layout()
        savefig(self.outputdir)

    def display_results(self):
        self._plot_matrix('corr', self.cols, self.corr, cmap='bwr')
        self._plot_matrix('pvals', self.cols, self.pvals, sum_rows=True)

        reorder_index = []
        with open(os.path.join(self.outputdir, 'ordered_pvals.csv'), 'w') as f:
            f.write('{},{}\n'.format('var', 'sum_of_pvals'))
            for i, c, p in sorted(zip(range(len(self.cols)), self.cols, self.pvals.sum(axis=0)),
                               key=lambda item: item[2]):
                f.write('{},{}\n'.format(c, p))
                reorder_index.append(i)

        ordered_corr = self.corr[reorder_index][:, reorder_index]
        ordered_pvals = self.pvals[reorder_index][:, reorder_index]
        ordered_cols = np.array(self.cols)[reorder_index]

        self._plot_matrix('ordered_corr', ordered_cols, ordered_corr, cmap='bwr')
        self._plot_matrix('ordered_pvals', ordered_cols, ordered_pvals, sum_rows=True)

