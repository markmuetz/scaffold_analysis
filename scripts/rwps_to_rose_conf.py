"""Script to extract RWPs from latest COSAR data (5/2/19) and save as rose-app-*.conf

Uses RWPs - remapped to have same labelling as in COSAR paper - and information from
high-resolution scaffold dumps (uses S0W0Forced) to work out pressures at each height.
Uses this info to interp pressures of RWPs to heights, and then writes out results
to OUTPUT_DIR as rose-app-*.conf files.
"""
import os
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import iris
import pandas as pd

import omnium as om
from omnium.consts import kappa, p_ref

BASE_DUMP_DIR = '/home/markmuetz/mirrors/archer/work/cylc-run/u-be530/share/data/history'
BASE_RWP_DIR = ('/home/markmuetz/mirrors/rdf/um10.9_runs/archive/u-au197/omnium_output/'
                'om_v0.11.1.0_cosar_v0.8.3.0_e889d0f4f8/P5Y_DP20')
OUTPUT_DIR = '/home/markmuetz/Dropbox/PhD/Results/rwp_rose_confs'


def get_pressures(filename):
    """Get z, pressure profile from a dump file by horizontal mean"""
    expt = os.path.basename(os.path.dirname(filename))
    cache_name = os.path.join(OUTPUT_DIR, expt + '_pressures.hdf')
    if not os.path.exists(cache_name):
        da = iris.load(filename)
        exnerp = om.utils.get_cube(da, 0, 255)
        z = exnerp.coord('atmosphere_hybrid_height_coordinate').points
        pdata = exnerp.data ** (1 / kappa) * p_ref
        p_profile = pdata.mean(axis=(1, 2))
        df = pd.DataFrame(data=p_profile[None, :], columns=z)
        df.to_hdf(cache_name, 'pressures')
    else:
        print('Using cached pressures')
        df = pd.read_hdf(cache_name)
        z = df.columns.values
        p_profile = df.values[0]
    return z, p_profile


def get_all_pressures():
    """Get pressures from all dump files in BASE_RWP_DIR"""
    dirs = sorted(glob(os.path.join(BASE_DUMP_DIR, 'S*')))
    zs = []
    p_profiles = []
    for d in dirs:
        filename = os.path.join(d, 'atmosa_da480.nc')
        if os.path.islink(filename):
            continue

        z, p_profile = get_pressures(filename)

        zs.append(z)
        p_profiles.append(p_profile)

    p_profiles = np.array(p_profiles)
    return zs, p_profiles


def plot_all_pressures(zs, p_profiles):
    """Plot all pressures vs z"""
    for z, p_profile in zip(zs, p_profiles):
        plt.plot(p_profile, z)

    plt.show()


def extract_rwps():
    """Extract RWPs from denormalized magnitudes and remapped labels"""
    # Using Remapped labels means thay are same as in COSAR paper.
    df_denorm_mag = pd.read_hdf(os.path.join(BASE_RWP_DIR, 'denorm_mag.hdf'))
    df_remapped_labels = pd.read_hdf(os.path.join(BASE_RWP_DIR, 'remapped_kmeans_labels.hdf'))
    labels = df_remapped_labels['nc-10_seed-391137'].values
    u = df_denorm_mag.values[:, :20]
    v = df_denorm_mag.values[:, 20:]

    rwps = []
    for i in range(10):
        rwp = np.concatenate((np.median(u[labels == i], axis=0),
                              np.median(v[labels == i], axis=0)))
        rwps.append(rwp)
    rwps = np.array(rwps)
    df_rwps = pd.DataFrame(data=rwps, index=np.arange(1, 11), columns=df_denorm_mag.columns)
    return df_rwps


def plot_rwps(df_rwps):
    """Plot RWPs as a sanity check"""
    p = list(range(50, 1050, 50))
    fig, axes = plt.subplots(2, 5, sharex=True)
    for i in range(10):
        u, v = df_rwps.values[i, :20], df_rwps.values[i, 20:]
        ax = axes.flatten()[i]
        ax.plot(u, p, 'b-')
        ax.plot(v, p, 'r-')
        ax.set_ylim((1000, 50))

    plt.show()


def output_hdf(df_rwps):
    df_rwps.to_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'pressure_rwps')


def load_rwps():
    return pd.read_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'))


def convert_rwp_pressures_to_heights():
    """Use info from dump for S0W0Forced to interp pressures to heights

    Decision, if pressure > max_pressure (pressure at 31.5 m), use 1000 hPa
    means lowest height in all RWPs is 0 m.
    """
    expt = 'S0W0Forced'
    cache_name = expt + '_pressures.hdf'

    df_p_profile = pd.read_hdf(os.path.join(OUTPUT_DIR, cache_name))
    z = df_p_profile.columns.values
    p_profile = df_p_profile.values[0] / 100

    df_rwps = pd.read_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'pressure_rwps')
    p = np.array(range(50, 1050, 50))
    print(p)

    # right=0 sets lowest height to 0.
    pressure_heights = np.interp(p, p_profile[::-1], z[::-1], right=0)
    df_heights = pd.DataFrame(data=pressure_heights[None, :],
                              columns=df_rwps.columns[:20])
    df_heights.to_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'pressure_heights')

    df_height_rwps_u = pd.DataFrame(data=df_rwps.values[:, :20],
                                    columns=pressure_heights,
                                    index=np.arange(1, 11))
    df_height_rwps_v = pd.DataFrame(data=df_rwps.values[:, 20:],
                                    columns=pressure_heights,
                                    index=np.arange(1, 11))

    df_height_rwps_u.to_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'height_rwps_u')
    df_height_rwps_v.to_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'height_rwps_v')


def output_rwp_rose_confs():
    """Output rose-app-RWP_C*.conf for each RWP"""
    height_rwps_u = pd.read_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'height_rwps_u')
    height_rwps_v = pd.read_hdf(os.path.join(OUTPUT_DIR, 'rwps.hdf'), 'height_rwps_v')
    z = height_rwps_u.columns.values

    for i in range(10):
        fmt_dict = {}
        fmt_dict['num_heights'] = len(z)
        fmt_dict['heights'] = ','.join(['{:.3f}'.format(zv) for zv in z[::-1]])
        fmt_dict['u_data'] = ','.join(['{:.3f}'.format(u) for u in height_rwps_u.values[i][::-1]])
        fmt_dict['v_data'] = ','.join(['{:.3f}'.format(v) for v in height_rwps_v.values[i][::-1]])

        fn = 'rose-app-RWP_C{}.conf'.format(i + 1)
        with open(os.path.join(OUTPUT_DIR, fn), 'w') as f:
            f.write('[namelist:recon_idealised]\n'
                    'num_uv_init_heights={num_heights}\n'
                    'u_init_data={u_data}\n'
                    'uv_init_height={heights}\n'
                    'v_init_data={v_data}\n'
                    '\n'
                    '[namelist:idealised]\n'
                    'num_uv_relax_heights={num_heights}\n'
                    'u_relax_data={u_data}\n'
                    'uv_relax_height={heights}\n'
                    'v_relax_data={v_data}\n'.format(**fmt_dict))


if __name__ == '__main__':
    try:
        p_profiles
        zs
    except NameError:
        zs, p_profiles = get_all_pressures()
    plot_all_pressures(zs, p_profiles)

    df_rwps = extract_rwps()
    plot_rwps(df_rwps)
    output_hdf(df_rwps)

    convert_rwp_pressures_to_heights()

    output_rwp_rose_confs()
