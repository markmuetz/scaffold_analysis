import os
from configparser import ConfigParser
import iris

BASEDIR = '/home/markmuetz/archer_mirror/nerc/um10.7_runs/postproc/u-ap347'

if __name__ == '__main__':
    for expt, fmt in [('S0', 'b-'), 
                      ('S4', 'g-'),
                      ('S0_1km', 'b--'),
                      ('S4_1km', 'g--')]:
        d = iris.load(os.path.join(BASEDIR, 'share/data/history/', expt, 'atmosa_da240.nc'))
        cp = ConfigParser()
        cp.add_section('namelist:idealise')

        th = d[71]
        mv = d[57]

        th_profile = th.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        mv_profile = mv.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        alt = th_profile.coord('altitude').points

        th_val = ['{:.2f}'.format(v) for v in th_profile.data]
        mv_val = ['{:.3e}'.format(v) for v in mv_profile.data]
        alt_val = ['{:.2f}'.format(v) for v in alt]

        print(','.join(th_val))
        print(','.join(mv_val))
        print(','.join(alt_val))

        cp.set('namelist:idealise', 'mv_init_data', ','.join(mv_val))
        cp.set('namelist:idealise', 'mv_init_height', ','.join(alt_val))

        cp.set('namelist:idealise', 'theta_init_data', ','.join(th_val))
        cp.set('namelist:idealise', 'theta_init_height', ','.join(alt_val))
        with open('rose-app-{}_init.conf'.format(expt), 'w') as f:
            cp.write(f)
