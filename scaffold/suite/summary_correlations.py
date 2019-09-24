import os
from logging import getLogger

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import iris

from omnium import Analyser
from omnium.utils import cd, get_cube_from_attr, cm_to_inch, latex_sigfig
from scaffold.expt_settings import cycle


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

NAME_MAP = {
    'LLS': 'LLS (s$^{-1}$)',
    'MLS': 'MLS (s$^{-1}$)',
    'LLWD': 'LLWD (m s$^{-1}$)',
    'MLWD': 'MLWD (m s$^{-1}$)',
    'mean_surf_wind': 'mean surf. wind (m s$^{-1}$)',
    'CAPE': 'CAPE (J kg$^{-1}$)',
    'all_lifetime': 'all cloud lifetime (min)',
    'simple_lifetime': 'simple cloud lifetime (min)',
    'complex_lifetime': 'complex cloud group lifetime (min)',
    'total_MF': 'total MF',
    'sigma': '$\sigma$',
    'cluster_index': 'cluster-index (km)',
    'suppr_index': 'suppr-index (km)',
    'suppr_radius': 'suppr-radius (km)',
    'cluster_radius': 'cluster-radius (km)',
    'mean_num_clds': 'mean num. clouds',
    'Tsurf': 'T. surf.',
}

CORR_NAME_MAP = {
    'LLS': 'LLS',
    'MLS': 'MLS',
    'LLWD': 'LLWD',
    'MLWD': 'MLWD',
    'mean_surf_wind': 'mean surf. wind',
    'CAPE': 'CAPE',
    'all_lifetime': 'all lifetime',
    'simple_lifetime': 'simple lifetime',
    'complex_lifetime': 'complex lifetime',
    'total_MF': 'total MF',
    'MF_per_cloud': 'MF per cloud',
    'sigma': '$\sigma$',
    'cluster_index': 'cluster-index',
    'suppr_index': 'suppr-index',
    'suppr_radius': 'suppr-radius',
    'cluster_radius': 'cluster-radius',
    'mean_num_clds': 'mean num. clouds',
    'Tsurf': 'T. surf.',
}


def plot_corr(df, c1, c2, ax=None):
    title = '{} vs {}'.format(c1, c2).replace(' ', '_')
    if not ax:
        plt.figure(title, figsize=cm_to_inch(16, 12))
        plt.clf()
        ax = plt.gca()

    d1, d2 = df[c1].values, df[c2].values
    if c1 == 'cluster_index':
        d1 = d1 / 1000 # Convert from m to km
    if c2 == 'cluster_index':
        d2 = d2 / 1000 # Convert from m to km
    # plt.title(title)


    markers = {'S': 'o', 'R': 'x'}
    for expt, v1, v2 in zip(df.expt, d1, d2):
        ax.scatter(v1, v2, color='k', marker=markers[expt[0]])
    lr = stats.linregress(d1, d2)
    x = np.array([d1.min(), d1.max()])
    y = lr.slope * x + lr.intercept

    ax.plot(x, y, 'k--', label='m={}\nc={}\nr$^2$={}\np={}'.format(latex_sigfig(lr.slope),
                                                                   latex_sigfig(lr.intercept),
                                                                   latex_sigfig(lr.rvalue ** 2),
                                                                   latex_sigfig(lr.pvalue)))
    ax.legend()
    ax.set_xlabel(NAME_MAP.get(c1, c1))
    ax.set_ylabel(NAME_MAP.get(c2, c2))
    # plt.tight_layout()

    return lr


def get_mass_flux(expt):
    mass_flux_combined_fn = '{}/atmos.mass_flux_combined.nc'.format(expt)
    logger.debug('loading mass_flux: {}', mass_flux_combined_fn)
    mfc = iris.load(mass_flux_combined_fn)
    mf = get_cube_from_attr(mfc, 'omnium_cube_id', 'mass_flux_z1_w1_qcl1')
    tmf = get_cube_from_attr(mfc, 'omnium_cube_id', 'total_mass_flux_z1_w1_qcl1')
    sigma = get_cube_from_attr(mfc, 'omnium_cube_id', 'sigma_z1_w1_qcl1')
    num_cld_snapshots = 10 * 48  # 10 days, every half hour.
    return {'mean_num_clds': mf.shape[0] / num_cld_snapshots, 'MF_per_cloud': mf.data.mean(),
            'total_MF': tmf.data.mean(), 'sigma': sigma.data.mean()}


def get_cloud_lifetimes(expt):
    cloud_tracking_fn = '{}/atmos.cloud_tracking_z1_t1.txt'.format(expt)
    logger.debug('loading cloud_tracking: {}', cloud_tracking_fn)
    with open(cloud_tracking_fn, 'r') as f:
        lines = f.readlines()

    simple_lifetime = float(lines[4].split(',')[-1].strip())
    complex_lifetime = float(lines[8].split(',')[-1].strip())
    all_lifetime = float(lines[10].split(',')[-1].strip())
    return {'simple_lifetime': simple_lifetime,
            'complex_lifetime': complex_lifetime,
            'all_lifetime': all_lifetime}


def savefig(cwd):
    fig_fn = '{}/summary_info_{}.png'.format(cwd, plt.gcf().get_label())
    logger.debug('saving fig to: {}'.format(fig_fn))
    plt.savefig(fig_fn)


class SummaryCorrelations(Analyser):
    analysis_name = 'summary_correlations'
    single_file = True

    # Completely ignore omnium file loading mechanism!
    input_dir = 'omnium_output/{version_dir}'
    input_filename = '{input_dir}/summary_correlations.dummy'
    output_dir = 'omnium_output/{version_dir}/suite_summary_correlations'
    output_filenames = ['{output_dir}/atmos.summary_correlations.dummy',
                        '{output_dir}/all_col_values.hdf']

    def load(self):
        self.basedir = os.path.dirname(self.task.filenames[0])
        self.outputdir = os.path.dirname(self.task.output_filenames[0])

        with cd(self.basedir):
            org_data_fn = ('suite_S0W0Forced_S4W0Forced_S0W5Forced_S4W5Forced_'
                           'RWP_C1_RWP_C2_RWP_C3_RWP_C4_RWP_C5_'
                           'RWP_C6_RWP_C7_RWP_C8_RWP_C9_RWP_C10'
                           '/org_data.csv')
            logger.debug('loading org_data: {}', org_data_fn)

            self.df_org = pd.read_csv(org_data_fn)
            self.df_tds = {}
            self.df_dyns = {}
            for expt in SORTED_EXPTS:
                si_fn = '{}/atmos.restart_dump_summary_info.hdf'.format(expt)
                logger.debug('loading summary_info: {}', si_fn)
                self.df_tds[expt] = pd.read_hdf(si_fn, 'thermodynamic_summary')
                self.df_dyns[expt] = pd.read_hdf(si_fn, 'dynamic_summary')

            self.cloud_lifetimes = {}
            self.mass_flux = {}
            for expt in SORTED_EXPTS:
                self.cloud_lifetimes[expt] = get_cloud_lifetimes(expt)
                self.mass_flux[expt] = get_mass_flux(expt)

    def run(self):
        for expt in SORTED_EXPTS:
            logger.debug(expt)
            df_td = self.df_tds[expt]
            # Skip first column: runid.
            for col in df_td.columns[1:]:
                logger.debug('  {}: {}, {}'.format(col, df_td[col].mean(), df_td[col].std()))
            df_dyn = self.df_dyns[expt]
            # Skip first column: runid.
            for col in df_dyn.columns[1:]:
                logger.debug('  {}: {}, {}'.format(col, df_dyn[col].mean(), df_dyn[col].std()))

        cols = (list(df_dyn.columns[1:]) + list(df_td.columns[1:]) +
                list(self.df_org.columns[2:]) +
                ['mean_num_clds', 'MF_per_cloud', 'total_MF', 'sigma'] +
                ['all_lifetime', 'simple_lifetime', 'complex_lifetime'])

        all_col_values = {'expt': SORTED_EXPTS}
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
                elif col in ['mean_num_clds', 'MF_per_cloud', 'total_MF', 'sigma']:
                    col_values.append(self.mass_flux[expt][col])
                elif col in ['all_lifetime', 'simple_lifetime', 'complex_lifetime']:
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
                logger.info('{} vs {}: r2-val={}, p-val={}', col1, col2, lr.rvalue**2, lr.pvalue)
                # axes[i, j].plot(d1, d2, 'kx')
        # savefig(cwd)

        self.cols = cols
        self.all_col_values = all_col_values
        self.corr = corr
        self.pvals = pvals
        self.df_all_col_values = pd.DataFrame(self.all_col_values)

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
        self.df_all_col_values.to_hdf(self.task.output_filenames[1], 'all_col_values')

    def _plot_matrix(self, name, cols, data, cmap=None, sum_rows=False):
        plt.figure(name, figsize=cm_to_inch(21, 21.0))
        plt.clf()
        if cmap:
            plt.imshow(data, cmap=cmap)
        else:
            plt.imshow(data)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                plt.annotate('{:.2f}'.format(data[i, j]), (i - 0.4, j + 0.1),
                             color='w', fontsize=8)

        if sum_rows:
            for j in range(data.shape[1]):
                plt.annotate('{:.2f}'.format(data[:, j].sum()), (data.shape[0] - 0.3, j + 0.1),
                             color='k', fontsize=8, annotation_clip=False)

        plt.xticks(range(len(cols)), [CORR_NAME_MAP.get(c, c) for c in cols], rotation=90)
        plt.yticks(range(len(cols)), [CORR_NAME_MAP.get(c, c) for c in cols])
        ax = plt.gca()
        xticklabels = ax.get_xticklabels()
        yticklabels = ax.get_yticklabels()

        def colour_ticks(ticklabels):
            for i, ticklab in enumerate(ticklabels):
                if i >= 17:
                    ticklab.set_color(cycle[2])
                elif i >= 9:
                    ticklab.set_color(cycle[0])
                elif i >= 3:
                    ticklab.set_color(cycle[1])

        colour_ticks(xticklabels)
        colour_ticks(yticklabels)

        plt.tight_layout()
        savefig(self.outputdir)

    def display_results(self):
        self._plot_corr_matrices()
        self._plot_indiv_corr()

    def _plot_corr_matrices(self):
        self._plot_matrix('corr', self.cols, self.corr, cmap='bwr')
        self._plot_matrix('pvals', self.cols, self.pvals, sum_rows=True)

        robust_corr = self.corr.copy()
        robust_corr[self.pvals > 0.01] = 0
        self._plot_matrix('robust_corr', self.cols, robust_corr, cmap='bwr')

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

        triu_pvals = self.pvals[np.triu(np.ones_like(self.pvals, dtype=bool), 1)]
        fig = plt.figure('pvals_hist')
        plt.hist(triu_pvals, bins=20, range=(0, 1))
        savefig(self.outputdir)

    def _plot_indiv_corr(self):
        corr_pairs = [
            ['LLS', 'cluster_index'],
            ['mean_surf_wind', 'CAPE'],
        ]

        for c1, c2 in corr_pairs:
            plot_corr(self.df_all_col_values, c1, c2)
            plt.tight_layout()
            savefig(self.outputdir)

        title = '{} vs {}'.format('cluster_index', 'lifetimes').replace(' ', '_')
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True,
                                            figsize=cm_to_inch(20, 10), num=title)
        plot_corr(self.df_all_col_values, 'cluster_index', 'all_lifetime', ax1)
        ax1.set_ylim((0, None))
        plot_corr(self.df_all_col_values, 'cluster_index', 'simple_lifetime', ax2)
        ax2.set_ylim((0, 70))
        plot_corr(self.df_all_col_values, 'cluster_index', 'complex_lifetime', ax3)
        ax3.set_ylim((0, None))
        plt.tight_layout()
        savefig(self.outputdir)

