""" Functions to analyze the sensor data. """

import numpy as np
import numpy.ma as ma


def estimate_g(sensor):
    """ Estimate g vector from acceleration data in iss system.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    if not "iss_ax" in sensor.data:
        raise KeyError(
            f"Sensor {sensor.name} does not have acceleration data in initial sensor system.")

    # flag times when foot is not on ground
    a = sensor.data["a"]
    delta = np.abs(a - 1)
    mask = delta > 0.01

    for d in ["x", "y", "z"]:
        sensor.data["iss_a" + d +
                    "_rest"] = ma.masked_array(sensor.data["iss_a" + d], mask=mask)
        # sensor.data["iss_a" + d + "_rest"][mask] = np.nan

    gx = np.average(sensor.data["iss_ax_rest"], )
    gy = np.average(sensor.data["iss_ay_rest"])
    gz = np.average(sensor.data["iss_az_rest"])

    sensor.data["iss_g"] = np.array([gx, gy, gz])


def remove_g(sensor):
    """ Remove g vector from acceleration data in iss system.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    estimate_g(sensor)
    for n, d in enumerate(["x", "y", "z"]):
        sensor.data["iss_a" + d + "_gr"] = sensor.data["iss_a" +
                                                       d] - sensor.data["iss_g"][n]
