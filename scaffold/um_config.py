"""UM run configuration read from an expt's rose-app-run.conf (replaces omnium's Expt)."""
from configparser import ConfigParser

MODEL_TYPES = {1: 'global_model',
               2: 'LAM_classic',
               3: 'LAM_EW_cyclic',
               4: 'LAM_bicyclic',
               5: 'SCM',
               6: 'SSFM'}


class UMConfig:
    def __init__(self, rose_app_run_conf_file):
        self.filename = rose_app_run_conf_file
        cp = ConfigParser()
        with open(rose_app_run_conf_file, 'r') as f:
            split_first_line = f.readline().split('=')
            assert split_first_line[0] == 'meta'
            # version looks like: um-atmos/vn11.0
            version = split_first_line[1].strip()
            self.um_version = tuple([int(v) for v in version.split('/')[-1][2:].split('.')])
            if self.um_version >= (10, 9):
                self._namelist_idealised = 'namelist:recon_idealised'
            else:
                self._namelist_idealised = 'namelist:idealise'
            second_line = f.readline()
            assert second_line.strip() == ''
            cp.read_file(f)
        self.um_config = cp

    @property
    def model_type(self):
        return MODEL_TYPES[self.um_config['namelist:model_domain'].getint('model_type')]

    @property
    def dx(self):
        return self.um_config[self._namelist_idealised].getfloat('delta_xi1')

    @property
    def dy(self):
        return self.um_config[self._namelist_idealised].getfloat('delta_xi2')

    @property
    def dt(self):
        return (self.um_config['namelist:nlstcgen'].getfloat('secs_per_periodim') /
                self.um_config['namelist:nlstcgen'].getfloat('steps_per_periodim'))

    @property
    def nx(self):
        return self.um_config['namelist:nlsizes'].getint('global_row_length')

    @property
    def ny(self):
        return self.um_config['namelist:nlsizes'].getint('global_rows')

    @property
    def lx(self):
        return self.dx * self.nx

    @property
    def ly(self):
        return self.dy * self.ny
