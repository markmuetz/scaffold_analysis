import os

try:
    import f90nml
except ImportError:
    pass
import numpy as np


class VertLev(object):
    def __init__(self, suite_dir):
        vertlevs_filename = os.path.join(suite_dir, 'app/um/file/rce_vertlevs.nml')

        self.vertlevs = f90nml.read(vertlevs_filename)['vertlevs']
        self.eta_theta = np.array(self.vertlevs['eta_theta'])
        self.eta_rho = np.array(self.vertlevs['eta_rho'])
        self.z_top = self.vertlevs['z_top_of_model']
        assert len(self.eta_theta) == len(self.eta_rho) + 1

        self.z_theta = self.eta_theta * self.z_top
        self.z_rho = self.eta_rho * self.z_top
        self.dz_theta = self.z_theta[1:] - self.z_theta[:-1]
