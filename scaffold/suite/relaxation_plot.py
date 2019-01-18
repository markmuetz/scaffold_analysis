import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser
from omnium.utils import get_cube


class RelaxationPlot(Analyser):
    analysis_name = 'relaxation_plot'
    multi_expt = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmos.456.pp2.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.relaxation_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        plt.clf()
        plt.xlabel('theta_inc (K day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            theta_inc = get_cube(self.expt_cubes[expt], 53, 181)

            z = theta_inc.coord('level_height').points
            colour = None

            for i in range(8):
                print(i)
                # 2880 steps_per_perdiodm: converts to K/day
                theta_inc_profile = theta_inc[i].data.mean(axis=(1, 2)) * 2880
                if i == 0:
                    plot = plt.plot(theta_inc_profile, z / 1000, label=expt)
                    colour = plot[0].get_color()
                else:
                    plt.plot(theta_inc_profile, z / 1000, color=colour)

        plt.legend()
        plt.savefig(self.file_path('theta_incs.png'))

        plt.clf()
        plt.xlabel('mv_inc (g kg$^{-1}$ day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            mv_inc = get_cube(self.expt_cubes[expt], 53, 182)
            z = mv_inc.coord('level_height').points
            colour = None

            for i in range(8):
                print(i)
                # 2880 steps_per_perdiodm: converts to g/kg/day
                mv_inc_profile = mv_inc[i].data.mean(axis=(1, 2)) * 2880 * 1000
                if i == 0:
                    plot = plt.plot(mv_inc_profile, z / 1000, label=expt)
                    colour = plot[0].get_color()
                else:
                    plt.plot(mv_inc_profile, z / 1000, color=colour)

        plt.legend()
        plt.savefig(self.file_path('mv_incs.png'))


