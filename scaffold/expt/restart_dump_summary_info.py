import os
from logging import getLogger

import iris
import pandas as pd

has_metpy = False
try:
    import metpy
    import metpy.calc as mpcalc
    from metpy.plots import SkewT
    from metpy.units import units
    has_metpy = True
except ImportError:
    pass

from omnium import Analyser
from omnium.consts import kappa, p_ref
from omnium.utils import get_cube

logger = getLogger('scaf.rdsa')
if not has_metpy:
    logger.warning('metpy not available')

class RestartDumpSummaryInfo(Analyser):
    """

    """
    analysis_name = 'restart_dump_summary_info'
    multi_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da???.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.restart_dump_summary_info.hdf']

    def load(self):
        self.dumps = {}
        for filename in self.task.filenames:
            basename = os.path.basename(filename)
            runid = int(basename[9:12])
            if self.settings.start_runid < runid <= self.settings.end_runid:
                logger.debug('adding runid: {}'.format(runid))
                dump = iris.load(filename)
                self.dumps[runid] = dump
            else:
                logger.debug('skipping runid: {}'.format(runid))

    def run(self):
        self.summary_info_data = []

        for runid, dump in self.dumps.items():
            self._calc_thermodynamic_summary_info(runid, dump)

        self.df_summary_info = pd.DataFrame(self.summary_info_data,
                                            columns=['runid', 'Tsurf',
                                                     'LCL', 'LFC', 'LNB',
                                                     'CAPE', 'CIN'])

    def _calc_thermodynamic_summary_info(self, runid, dump):
        theta = get_cube(dump, 0, 4)
        exnerp = get_cube(dump, 0, 255)
        qv = get_cube(dump, 0, 10)

        qvdata = qv.data

        Tdata = theta.data * exnerp.data
        pdata = exnerp.data ** (1 / kappa) * p_ref

        p = pdata * units('Pa')
        qv = qvdata * units('kg/kg')
        T = Tdata * units('K')
        Td = mpcalc.dewpoint_from_specific_humidity(qv, T, p)

        p_profile = p.mean(axis=(1, 2))
        T_profile = T.mean(axis=(1, 2))
        Td_profile = Td.mean(axis=(1, 2))

        cape, cin = mpcalc.surface_based_cape_cin(p_profile, T_profile, Td_profile)
        lcl = mpcalc.lcl(p_profile[0], T_profile[0], Td_profile[0])
        lfc = mpcalc.lfc(p_profile, T_profile, Td_profile)
        lnb = mpcalc.el(p_profile, T_profile, Td_profile)

        # logger.info(lfc)
        # logger.info(mpcalc.pressure_to_height_std(lfc[0]))
        self.summary_info_data.append([runid,
                                       T_profile[0].magnitude,
                                       lcl[0].magnitude,
                                       lfc[0].magnitude,
                                       lnb[0].magnitude,
                                       cape.magnitude,
                                       cin.magnitude])



    def save(self, state, suite):
        self.df_summary_info.to_hdf(self.task.output_filenames[0], 'thermodynamic_summary')

