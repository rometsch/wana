""" Integration of acceleration to estimate velocities and positions. """
import numpy as np
from wana.analysis import calculate_norm

def estimate_velocities(sensor, frame, perstep=False):
    """ Use accelerations with g removed to estimate velocities in the iss.

    iss = initial sensor system
    lab = laboratory frame
    Integrate the accelerations.

    Adds arrays containing the velocity: "{iss,lab}_v{x,y,z}".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    frame: str
        Lab or iss frame.
    perstep: bool
        Reset to zero before each new step.
    """
    dt = sensor.data["dt"]
    for d in ["x", "y", "z"]:
        acc_name = frame + "_a" + d + "_gr"
        a = sensor.data[acc_name]

        # second order integration
        aux = np.zeros(len(dt))
        aux[1:] = (a[1:] + a[:-1])/2
        aux[0] = a[0]
        aux *= dt

        v = np.cumsum(aux)

        if perstep:
            index_step_start = sensor.data["interval_steps"][:, 0]
            for i in index_step_start:
                v[i:] -= v[i]

        varname = frame + "_v" + d
        if perstep:
            varname += "_step"
        sensor.data[varname] = v
        sensor.units[varname] = "m/s"

    varpattern = frame + "_v{}"
    if perstep:
        varpattern += "_step"
    calculate_norm(sensor, varpattern, unit="m/s")


def estimate_positions(sensor, frame, varpattern, perstep=False):
    """ Use velocities to estimate positions in the iss.

    iss = initial sensor system
    Integrate the velocities.

    Adds arrays containing the position: "{iss,lab}_{x,y,z}".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    frame: str
        Lab or iss frame.
    perstep: bool
        Reset to zero before each new step.
    """
    dt = sensor.data["dt"]
    for d in ["x", "y", "z"]:
        vel_name = frame + "_" + varpattern.format(d)
        v = sensor.data[vel_name]
        # x = np.cumsum(v*dt)

        # second order integration
        aux = np.zeros(len(dt))
        aux[1:] = (v[1:] + v[:-1])/2
        aux[0] = v[0]
        aux *= dt

        x = np.cumsum(aux)

        if perstep:
            index_step_start = sensor.data["interval_steps"][:, 0]
            for i in index_step_start:
                x[i:] -= x[i]

        varname = frame + "_" + varpattern[1:].format(d)
        if perstep:
            varname += "_step"
        sensor.data[varname] = x
        sensor.units[varname] = "m"

    varpattern = frame + "_" + varpattern[1:]
    calculate_norm(sensor, varpattern, unit="m")


def linear_step_correction(sensor, varname, unit=None):
    """ Use values at beginning and end of step to correct for drift.

    Use a linear interpolation to compute the correction.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    varname: str
        Variable name.
    unit: str
        Physical unit of the variable.
    """
    step_intervals = sensor.data["interval_steps"]

    y = sensor.data[varname]

    yc = np.zeros(len(y))

    for I in step_intervals:
        left = I[0]
        right = I[1]

        dy = y[right] - y[left]
        dx = np.linspace(0, 1, num=right-left)

        correction = - y[left] - dx*dy

        yc[left:right] = y[left:right] + correction

    sensor.data[varname + "_lc"] = yc

    if unit is not None:
        sensor.units[varname + "_lc"] = unit
