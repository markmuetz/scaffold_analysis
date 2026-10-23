"""metpy calls whose signatures changed between metpy 0.12 (omnium env) and 1.x (remake env)."""
import metpy
import metpy.calc as mpcalc


def dewpoint_from_specific_humidity(specific_humidity, temperature, pressure):
    """Dewpoint, whatever the installed metpy's argument order.

    metpy < 1.0: (specific_humidity, temperature, pressure); 1.0-1.5: (pressure, temperature,
    specific_humidity); >= 1.6: (pressure, specific_humidity) -- temperature is not needed.
    """
    version = tuple(int(v) for v in metpy.__version__.split('.')[:2])
    if version < (1, 0):
        return mpcalc.dewpoint_from_specific_humidity(specific_humidity, temperature, pressure)
    if version < (1, 6):
        return mpcalc.dewpoint_from_specific_humidity(pressure, temperature, specific_humidity)
    return mpcalc.dewpoint_from_specific_humidity(pressure, specific_humidity)
