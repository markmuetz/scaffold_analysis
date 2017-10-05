import os
import pylab as plt
from configparser import ConfigParser

from settings import BASEDIR, OUTDIR

def plot_from_file(expt, filename, fmt):
    cp = ConfigParser()
    if expt == 'init':
        with open(filename, 'r') as f:
            f.readline()
            f.readline()
            cp.read_file(f)
    else:
        cp.read(filename)

    ideal = cp['namelist:idealise']

    mv_init_height = []
    mv_init_data = []
    th_init_height = []
    th_init_data = []

    for mv, mv_h in zip(ideal['mv_init_data'].split(','), ideal['mv_init_height'].split(',')):
        mv_init_data.append(float(unicode.replace(mv, '\n=', '')))
        mv_init_height.append(float(unicode.replace(mv_h, '\n=', '')))
        
        
    for th, th_h in zip(ideal['theta_init_data'].split(','), ideal['theta_init_height'].split(',')):
        th_init_data.append(float(unicode.replace(th, '\n=', '')))
        th_init_height.append(float(unicode.replace(th_h, '\n=', '')))
        
    plt.figure('mv')
    plt.plot(mv_init_data, mv_init_height, fmt, label=expt)
    plt.figure('th')
    plt.plot(th_init_data, th_init_height, fmt, label=expt)


def main():
    filename = os.path.join(BASEDIR, 'app/um/rose-app.conf')
    plot_from_file('init', filename, 'k:')

    for expt, fmt in [('S0', 'b-'), 
                      ('S4', 'g-'),
                      ('S0_1km', 'b--'),
                      ('S4_1km', 'g--')]:
        filename = os.path.join(OUTDIR, 'rose-app-{}_init.conf'.format(expt))
        plot_from_file(expt, filename, fmt)

    plt.figure('mv')
    plt.legend()
    plt.figure('th')
    plt.legend()

    plt.show()


if __name__ == '__main__':
    main()
