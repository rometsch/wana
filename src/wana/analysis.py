""" Functions to analyze the sensor data. """

import numpy as np
import numpy.ma as ma

import wana.event_detection as ed


def estimate_g(sensor, order=1):
    """ Estimate g vector from acceleration data in iss system.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    order: int
        Order of the polynomial to fit.
    """
    mask_a_resting = sensor.data["mask_a_resting"]
    for v in ["iss_ax", "iss_ay", "iss_az"]:
        ed.flag_resting(sensor, v, prior_mask=mask_a_resting,
                        N_kernel=1, delta_threshold=0.1,
                        order=order)
        extract_resting_values(sensor, v)
    estimate_g_poly(sensor, order)
    set_g_units(sensor)


def extract_resting_values(sensor, varname):
    """ Extract values for the times when sensor is at rest.

    iss = initial sensor system

    Adds masked arrays containing the accelerations when the sensor is at rest: "{varname}_resting"
    and a time array at corresponding indices: "time_{varname}_resting".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    time = sensor.data["time"]
    mask_resting = np.logical_not(sensor.data[f"mask_{varname}_resting"])

    sensor.data[f"time_{varname}_resting"] = ma.masked_array(
        time, mask=mask_resting)
    sensor.units[f"time_{varname}_resting"] = "s"

    values_masked = ma.masked_array(sensor.data[varname], mask=mask_resting)

    sensor.data[varname + "_resting"] = values_masked
    sensor.units[varname + "_resting"] = "m/s2"


def set_g_units(sensor):
    """ Set units for the g vector data.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    sensor.units["iss_g"] = "m/s2"
    sensor.units["iss_gx"] = "m/s2"
    sensor.units["iss_gy"] = "m/s2"
    sensor.units["iss_gz"] = "m/s2"


def estimate_g_average(sensor):
    """ Estimate g vector by averaging accelerations.

    iss = initial sensor system

    Adds a len 3 array for the vector of g (iss_g) and array for the
    components of g at every point in time (iss_g{x,y,z}).

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    gx = np.average(sensor.data["iss_ax_resting"])
    gy = np.average(sensor.data["iss_ay_resting"])
    gz = np.average(sensor.data["iss_az_resting"])

    g = np.sqrt(gx**2 + gy**2 + gz**2)

    sensor.data["iss_g"] = np.array([gx, gy, gz])  # /g*g_constant
    sensor.units["iss_g"] = "m/s2"

    # add components
    N = len(sensor.data["counter"])
    sensor.data["iss_gx"] = np.ones(N)*gx
    sensor.data["iss_gy"] = np.ones(N)*gy
    sensor.data["iss_gz"] = np.ones(N)*gz


def estimate_g_poly(sensor, order, mode="single"):
    """ Estimate g vector by fitting a polynomial to accelerations.

    iss = initial sensor system


    Adds arrays for the components of g at every point in time (iss_g{x,y,z}).

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    order: int
        Order of the polynomial.
    """
    time = sensor.data["time"]
    for d in ["x", "y", "z"]:
        varname = f"iss_a{d}"
        if mode == "single":
            mask = sensor.data[f"mask_{varname}_resting"]
        elif mode == "combined":
            mask = sensor.data[f"mask_a_resting"]
        else:
            raise ValueError("Mode", mode, "is not supported!")
        time_resting = time[mask]
        values_resting = sensor.data[varname][mask]

        p = np.polyfit(time_resting, values_resting, order)
        f = np.poly1d(p)
        g = f(time)
        sensor.data[f"iss_g{d}"] = g


def remove_g(sensor):
    """ Remove g vector from acceleration data in iss system.

    iss = initial sensor system

    Adds arrays containing the accelerations with g removed: "iss_a{x,y,z}_gr".

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    for n, d in enumerate(["x", "y", "z"]):
        g_n = sensor.data["iss_g"+d]
        a = sensor.data["iss_a"+d]

        varname = "iss_a" + d + "_gr"
        sensor.data[varname] = a - g_n
        sensor.units[varname] = "m/s2"

    calculate_norm(sensor, "iss_a{}_gr", unit="m")


def zero_acceleration_resting(sensor, frame):
    """ Set acceleration with g removed to zero in resting period.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    frame: str
        Coordinate system.
    """
    mask_resting = sensor.data["mask_a_resting"]
    for n, d in enumerate(["x", "y", "z"]):
        varname = frame + "_a" + d + "_gr"
        sensor.data[varname][mask_resting] = 0

    calculate_norm(sensor, frame+"_a{}_gr", unit="m")


def calculate_norm(sensor, varpattern, unit):
    """ Calculate absolute magnitude of vector.

    Add an array containing the values to data.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    varpattern: str
        Variable name.
    unit: str
        Physical unit of the variable.
    """
    x = sensor.data[varpattern.format("x")]
    y = sensor.data[varpattern.format("y")]
    z = sensor.data[varpattern.format("z")]

    l = np.sqrt(x**2 + y**2 + z**2)

    varname = varpattern.format("")
    sensor.data[varname] = x
    if unit is not None:
        sensor.units[varname] = unit


def calculate_velocity_direction_angles(sensor, varpattern):
    """ Calculate the direction of the velocity vector.

    Add an array containing the values to data.

    The varpattern includes a {} to be filled with the directions:
    e.g. lab_a{}_lc

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    varpattern: str
        Pattern of the variable.
    """
    vec = np.array([
        sensor.data[varpattern.format("x")],
        sensor.data[varpattern.format("y")],
        sensor.data[varpattern.format("z")]
    ])
    length = np.linalg.norm(vec, axis=0)

    angles = np.arccos(vec/length)

    sensor.data[varpattern.format("x")+"_ori"] = angles[0]
    sensor.data[varpattern.format("y")+"_ori"] = angles[1]
    sensor.data[varpattern.format("z")+"_ori"] = angles[2]
