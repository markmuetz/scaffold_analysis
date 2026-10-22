"""Loading/saving helpers that replace the omnium Analyser plumbing."""
import os
from collections import OrderedDict
from logging import getLogger

import iris

logger = getLogger('scaffold.io')


def load_cubes(filenames):
    """Load one file, or several files concatenated (omnium's single_file/multi_file loading)."""
    if isinstance(filenames, (str, os.PathLike)):
        return iris.load(str(filenames))
    cubes = iris.load([str(fn) for fn in filenames])
    return iris.cube.CubeList.concatenate(cubes)


def load_expt_cubes(expt_filenames):
    """Load one file per expt (omnium's multi_expt loading): {expt: cubes}."""
    return OrderedDict((expt, iris.load(str(fn))) for expt, fn in expt_filenames.items())


def save_cubes(results, filename, analysis_name, zlib=True, complevel=None):
    """Save {cube_id: cube} to netCDF.

    Each cube gets the `omnium_cube_id` and `omnium_process` attributes (downstream analyses
    select cubes by `omnium_cube_id`), but no version/git attributes.
    """
    for cube_id, cube in results.items():
        cube.attributes['omnium_cube_id'] = cube_id
        cube.attributes['omnium_process'] = analysis_name
    cubelist = iris.cube.CubeList(list(results.values()))
    if not len(cubelist):
        logger.warning('No results to save')
        return
    save_kwargs = {'zlib': zlib}
    if complevel is not None:
        save_kwargs['complevel'] = complevel
    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    if os.path.islink(filename):
        os.remove(filename)
    iris.save(cubelist, str(filename), **save_kwargs)


def side_path(output_filename, analysis_name, name):
    """Path for an extra output next to a task's main output (omnium's Analyser.file_path)."""
    dirname = os.path.dirname(str(output_filename))
    filename = os.path.join(dirname, analysis_name + '_' + name)
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    return filename
