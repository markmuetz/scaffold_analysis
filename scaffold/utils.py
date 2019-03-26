import numpy as np


def cm_to_inch(*vals):
    return [v / 2.54 for v in vals]


def find_intersections(c1, c2):
    signs = np.ones_like(c1)
    diff = c1 - c2
    signs[diff < 0] = -1
    signs[diff == 0] = 0
    indices = np.where((np.abs(np.roll(signs, -1) - signs) == 2) | (signs == 0))[0]
    weights = np.array([1 / ((np.abs(diff[(i + 1) % len(diff)] / diff[i])) + 1) for i in indices])
    return indices, weights


def interp_vert_rho2w(vertlevs, w_slice, rho_slice, time_index, height_level_index,
                      level_number):
    """Interpolate vertically from rho to w. Perform checks to make sure interp. is good. """

    # Work out in 2 ways and check equal because paranoid.
    # AKA sanity checks.
    w_height = w_slice.attributes['heights'][height_level_index]
    w_height2 = vertlevs.z_theta[level_number]
    # N.B. for every 1 w_height, there are 2 rho_heights. Index appropriately.
    rho_height_lower = rho_slice.attributes['heights'][2 * height_level_index]
    rho_height_lower2 = vertlevs.z_rho[level_number - 1]
    rho_height_upper = rho_slice.attributes['heights'][2 * height_level_index + 1]
    rho_height_upper2 = vertlevs.z_rho[level_number]
    # Calc scaling for linear interp.
    alpha = (w_height - rho_height_lower) / (rho_height_upper - rho_height_lower)

    # Paranoia. Well justified it turns out. Saved me from doing wrong analysis.
    assert w_height == w_height2
    assert rho_height_lower == rho_height_lower2
    assert rho_height_upper == rho_height_upper2
    assert 0 <= alpha <= 1

    # Interp rho onto w grid.
    rho_ss_lower = rho_slice[time_index, height_level_index].data
    rho_ss_upper = rho_slice[time_index, height_level_index + 1].data

    rho_ss_interp = (1 - alpha) * rho_ss_lower + alpha * rho_ss_upper
    return rho_ss_interp
