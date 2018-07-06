import matplotlib
matplotlib.use('Agg')
import pylab as plt
from matplotlib.colors import LogNorm

from omnium.analyser import Analyser
from omnium.utils import get_cube


class PrecipPlot(Analyser):
    """Pick out precip timesteps and plot."""
    analysis_name = 'precip_plot'
    multi_expt = True
    input_dir = 'work/20000101T0000Z/{expt}_atmos'
    input_filename_glob = '{input_dir}/atmos.pp3.nc'
    output_dir = 'omnium_output/{version_dir}/suite'
    output_filenames = ['{output_dir}/atmos.precip_plot.dummy']

    expts_to_plot = None

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self._plot()

    def _plot(self):
        if not self.expts_to_plot:
            self.expts_to_plot = self.task.expts

        precips = {}
        for expt in self.task.expts:
            cubes = self.expt_cubes[expt]
            precip = get_cube(cubes, 4, 203)
            precips[expt] = precip

        max_precips = ['time_index, expt, max precip [mm/hr]']
        for i in range(precip.shape[0] - 100, precip.shape[0]):
            fig, axes = plt.subplots(1, len(self.expts_to_plot))

            precip_max = 0
            for expt in self.expts_to_plot:
                precip = precips[expt]
                precip_max = max(precip[i].data.max(), precip_max)

            for ax, expt in zip(axes, self.expts_to_plot):
                ax.set_title(expt)
                if expt == self.task.expts[0]:
                    ax.set_ylabel('y (km)')
                else:
                    ax.get_yaxis().set_visible(False)
                ax.set_xlabel('x (km)')

                precip = precips[expt]
                # precip in kg m-2 s-1, want mm hr-1:
                # /rho_water: m s-1
                # *1000: mm s-1
                # *3600: mm hr-1
                # N.B. rho_water = 1000 kg m-3.
                #import ipdb; ipdb.set_trace()
                precip_data = precip[i].data * 3600

                precip_min = 1e-4
                precip_data[precip_data < precip_min] = 0
                im = ax.imshow(precip_data, origin='lower',
                               interpolation='nearest', extent=[0, 256, 0, 256],
                               #vmin=0, vmax=precip_max * 3600)
                               norm=LogNorm(vmin=precip_min, vmax=precip_max * 3600))

            for expt in self.task.expts:
                precip = precips[expt]
                precip_data = precip[i].data * 3600
                max_precips.append('{},{},{}'.format(i, expt, precip_data.max()))

            plt.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.89, 0.27, 0.02, 0.46])
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label('rainfall (mm hr$^{-1}$)', rotation=270, labelpad=20)

            plt.savefig(self.file_path('time_index{}.png'.format(i)))
            plt.close('all')
        self.save_text('max_precip.csv', '\n'.join(max_precips) + '\n')
