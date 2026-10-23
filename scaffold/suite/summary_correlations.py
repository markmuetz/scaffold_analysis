import os
from logging import getLogger

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import iris

from omnium import Analyser
from omnium.utils import cd, get_cube_from_attr, cm_to_inch, latex_sigfig
from scaffold.expt_settings import cycle, EXPT_DETAILS
from scaffold.paper_style import use_paper_style, panel_labels, save_paper_fig, FULL_WIDTH_CM, SINGLE_COL_CM


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
    'LLWD': 'LLS (m s$^{-1}$)',
    'MLWD': 'MLS (m s$^{-1}$)',
    'mean_surf_wind': 'surf. wind (m s$^{-1}$)',
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
    'Tsurf': '$T_\\mathrm{s}$ (K)',
}

CORR_NAME_MAP = {
    'LLS': 'LLS',
    'MLS': 'MLS',
    'LLWD': 'LLS',
    'MLWD': 'MLS',
    'mean_surf_wind': 'surf. wind',
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
    'Tsurf': '$T_\\mathrm{s}$',
}

# Groups of variables in the correlation matrix, in column order.
CORR_GROUPS = [('Driving', 3), ('Thermodynamic', 6), ('Cloud field', 8), ('Lifetime', 3)]

# RWPs labelled directly in the scatter plots.
LABELLED_RWPS = ['RWP_C1', 'RWP_C2', 'RWP_C6', 'RWP_C10']


def _fmt_r(r):
    """Correlation to 2 d.p. without the leading zero, with a true minus sign: -0.94 -> −.94."""
    text = '{:.2f}'.format(abs(r))
    if text.startswith('0'):
        text = text[1:]
    return ('\u2212' if r < 0 else '') + text


def _fmt_num(v):
    return latex_sigfig(v).replace('-', '\u2212')


def _fmt_p(p):
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '{}$\\times 10^{{{}}}$'.format(mantissa, int(exponent)).replace('-', '\u2212')


def _empty_corner(x, y):
    """Pick the axes corner with the fewest points, for the regression annotation."""
    xn = (x - x.min()) / (np.ptp(x) or 1)
    yn = (y - y.min()) / (np.ptp(y) or 1)
    corners = {'upper left': ((xn < 0.5) & (yn > 0.5)).sum(),
               'upper right': ((xn > 0.5) & (yn > 0.5)).sum(),
               'lower right': ((xn > 0.5) & (yn < 0.5)).sum(),
               'lower left': ((xn < 0.5) & (yn < 0.5)).sum()}
    return min(corners, key=corners.get)


def plot_corr(df, c1, c2, ax=None, legend=True, direct_labels=True, annotate_xy=None,
              legend_loc='best', label_offsets=None):
    title = '{} vs {}'.format(c1, c2).replace(' ', '_')
    if not ax:
        plt.figure(title, figsize=cm_to_inch(SINGLE_COL_CM, 6.5))
        plt.clf()
        ax = plt.gca()

    d1, d2 = df[c1].values, df[c2].values
    if c1 == 'cluster_index':
        d1 = d1 / 1000 # Convert from m to km
    if c2 == 'cluster_index':
        d2 = d2 / 1000 # Convert from m to km

    # Circles = simple profiles (experiment colours), crosses = RWPs (black).
    for expt, v1, v2 in zip(df.expt, d1, d2):
        if expt[0] == 'S':
            ax.scatter(v1, v2, color=EXPT_DETAILS[expt][1], marker='o', s=18, zorder=3,
                       edgecolor='none')
        else:
            ax.scatter(v1, v2, color='k', marker='x', s=18, linewidth=0.8, zorder=3)
        if direct_labels and expt in LABELLED_RWPS:
            offset = (label_offsets or {}).get(expt, (3, 2))
            ax.annotate(EXPT_DETAILS[expt][0], (v1, v2), xytext=offset,
                        textcoords='offset points', fontsize=6.5,
                        ha='center' if offset[1] < 0 else ('left' if offset[0] >= 0 else 'right'))
    lr = stats.linregress(d1, d2)
    x = np.array([d1.min(), d1.max()])
    y = lr.slope * x + lr.intercept
    ax.plot(x, y, color='0.4', linestyle='--', linewidth=1, zorder=2)

    sign = '+' if lr.intercept >= 0 else '\u2212'
    text = 'y = {}x {} {}\nr$^2$ = {}\np = {}'.format(_fmt_num(lr.slope), sign,
                                                     _fmt_num(abs(lr.intercept)),
                                                     _fmt_num(lr.rvalue ** 2),
                                                     _fmt_p(lr.pvalue))
    corner = annotate_xy or _empty_corner(d1, d2)
    xa = 0.04 if 'left' in corner else 0.96
    # Leave room for a panel label at the top-left.
    ya = 0.86 if corner == 'upper left' else (0.96 if 'upper' in corner else 0.04)
    ax.text(xa, ya, text, transform=ax.transAxes, fontsize=7,
            ha='left' if 'left' in corner else 'right',
            va='top' if 'upper' in corner else 'bottom')

    if legend:
        handles = [Line2D([], [], marker='o', color='0.3', linestyle='none', markersize=4,
                          label='simple profiles'),
                   Line2D([], [], marker='x', color='k', linestyle='none', markersize=4,
                          label='RWPs')]
        ax.legend(handles=handles, loc=legend_loc, handletextpad=0.2)
    ax.set_xlabel(NAME_MAP.get(c1, c1))
    ax.set_ylabel(NAME_MAP.get(c2, c2))

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
    fig_fn = '{}/summary_info_{}.pdf'.format(cwd, plt.gcf().get_label())
    logger.debug('saving fig to: {}'.format(fig_fn))
    save_paper_fig(fig_fn)


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
        print('oh hai!!!')
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

        # Copy self.cols
        if 'LLS' in cols:
            cols.remove('LLS')
        if 'MLS' in cols:
            cols.remove('MLS')

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
        use_paper_style()
        self._plot_corr_matrices()
        self._plot_indiv_corr()

    def _plot_corr_matrices(self):
        self._plot_matrix('corr', self.cols, self.corr, cmap='bwr')
        self._plot_matrix('pvals', self.cols, self.pvals, sum_rows=True)

        self._plot_corr_triangle(self.cols, self.corr, self.pvals)

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

    def _plot_corr_triangle(self, cols, corr, pvals):
        """Paper Fig. 9: lower-triangle correlation matrix of the robust (p < 0.01) correlations.

        Cells passing a Bonferroni correction are outlined and their values are in bold.
        """
        n = len(cols)
        npairs = n * (n - 1) // 2
        lower = np.tril(np.ones((n, n), dtype=bool), -1)
        robust = lower & (pvals < 0.01)
        bonf_thresh = 0.05 / npairs
        bonf = lower & (pvals < bonf_thresh)
        logger.info('Bonferroni: {} pairs, {} pass p < {}', npairs, bonf.sum(), bonf_thresh)

        fig, ax = plt.subplots(num='robust_corr', figsize=cm_to_inch(FULL_WIDTH_CM, 15))
        ax.imshow(np.where(lower, 1., np.nan), cmap=ListedColormap(['#ececec']), vmin=0, vmax=1)
        im = ax.imshow(np.ma.masked_where(~robust, corr), cmap='RdBu_r', vmin=-1, vmax=1)
        for i, j in zip(*np.where(robust)):
            ax.text(j, i, _fmt_r(corr[i, j]), ha='center', va='center', fontsize=6,
                    color='white' if abs(corr[i, j]) > 0.6 else 'black',
                    fontweight='bold' if bonf[i, j] else 'normal')
            if bonf[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec='#1b7837',
                                           lw=1.5))

        group_edges = np.cumsum([g[1] for g in CORR_GROUPS])
        assert group_edges[-1] == n, 'groups do not match the columns'
        for e in group_edges[:-1] - 0.5:
            ax.plot([-0.5, e], [e, e], 'k', lw=0.6)
            ax.plot([e, e], [e, n - 0.5], 'k', lw=0.6)

        labels = [CORR_NAME_MAP.get(c, c) for c in cols]
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=7)
        ax.tick_params(length=0)
        # The first row and last column of the lower triangle are empty.
        ax.set_xlim(-0.5, n - 1.5)
        ax.set_ylim(n - 0.5, 0.5)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Group bracket labels on the left and bottom edges.
        starts = np.concatenate([[0], group_edges[:-1]])
        trans_y = blended_transform_factory(ax.transAxes, ax.transData)
        trans_x = blended_transform_factory(ax.transData, ax.transAxes)
        for (name, _), start, end in zip(CORR_GROUPS, starts, group_edges):
            # Rows shown are 1..n-1, columns 0..n-2.
            r0, r1 = max(start, 1) - 0.4, end - 0.6
            if r1 > r0:
                ax.plot([-0.30, -0.30], [r0, r1], color='k', lw=0.6, transform=trans_y,
                        clip_on=False)
                ax.text(-0.32, (r0 + r1) / 2, name, transform=trans_y, rotation=90,
                        ha='right', va='center', fontsize=7)
            c0, c1 = start - 0.4, min(end, n - 1) - 0.6
            if c1 > c0:
                ax.plot([c0, c1], [-0.24, -0.24], color='k', lw=0.6, transform=trans_x,
                        clip_on=False)
                ax.text((c0 + c1) / 2, -0.26, name, transform=trans_x, ha='center', va='top',
                        fontsize=7)

        fig.colorbar(im, ax=ax, shrink=0.6, label='$r$')
        savefig(self.outputdir)

    def _plot_indiv_corr(self):
        corr_pairs = [
            ['LLWD', 'cluster_index'],
            ['mean_surf_wind', 'CAPE'],
        ]

        annotate_corners = {('LLWD', 'cluster_index'): 'lower right'}
        for c1, c2 in corr_pairs:
            plot_corr(self.df_all_col_values, c1, c2,
                      annotate_xy=annotate_corners.get((c1, c2)))
            savefig(self.outputdir)

        title = '{} vs {}'.format('cluster_index', 'lifetimes').replace(' ', '_')
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True,
                                            figsize=cm_to_inch(FULL_WIDTH_CM, 6.5), num=title)
        # In (a) and (c) C1 sits among other crosses: put its label just below the point.
        offsets = {'RWP_C1': (0, -10)}
        plot_corr(self.df_all_col_values, 'cluster_index', 'all_lifetime', ax1, legend=False,
                  label_offsets=offsets)
        ax1.set_ylim((0, None))
        # One marker legend for the figure, in the empty lower left of panel (b).
        plot_corr(self.df_all_col_values, 'cluster_index', 'simple_lifetime', ax2,
                  legend_loc='lower left')
        ax2.set_ylim((0, 70))
        plot_corr(self.df_all_col_values, 'cluster_index', 'complex_lifetime', ax3, legend=False,
                  label_offsets=offsets)
        ax3.set_ylim((0, None))
        for ax in (ax1, ax2, ax3):
            # The text discusses the cluster-index threshold of 75 km.
            ax.axvline(75, color='0.8', lw=0.6, zorder=1)
        panel_labels((ax1, ax2, ax3))
        fig.tight_layout()
        savefig(self.outputdir)
