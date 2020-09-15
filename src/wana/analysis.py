""" Functions to analyze the sensor data. """

import numpy as np
import numpy.ma as ma

# the numerical constant of earth's acceleration
g_constant = 9.81


def flag_resting(sensor):
    """ Get a mask for when the sensor is resting on ground.

    Compare the total acceleration to the value of g = 9.81 m/s2

    Adds a bool array "mask_resting.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    a = sensor.data["a"]
    delta = np.abs(a - g_constant)
    mask = delta < 0.05*g_constant

    N_kernel = 10
    kernel = np.ones(N_kernel)/N_kernel
    conv = np.convolve(mask, kernel, mode="same")
    sensor.data["mask_resting"] = conv > 0.5
        


def estimate_g(sensor):
    """ Estimate g vector from acceleration data in iss system.

    iss = initial sensor system

    Adds masked arrays containing the accelerations when the sensor is at rest: "iss_a{x,y,z}_rest".

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    flag_resting(sensor)
    mask_moving = np.logical_not(sensor.data["mask_resting"])
    for d in ["x", "y", "z"]:
        a_masked = ma.masked_array(sensor.data["iss_a" + d], mask=mask_moving)
        
        varname = "iss_a" + d + "_rest"
        sensor.data[varname] = a_masked
        sensor.units[varname] = "m/s2"

    gx = np.average(sensor.data["iss_ax_rest"])
    gy = np.average(sensor.data["iss_ay_rest"])
    gz = np.average(sensor.data["iss_az_rest"])

    g = np.sqrt(gx**2 + gy**2 + gz**2)
    
    sensor.data["iss_g"] = np.array([gx, gy, gz])/g*g_constant
    sensor.units["iss_g"] = "m/s2"
    
    print(f"g = ({gx}, {gy}, {gz}), |g| = {g}")


def remove_g(sensor):
    """ Remove g vector from acceleration data in iss system.

    iss = initial sensor system

    Adds arrays containing the accelerations with g removed: "iss_a{x,y,z}_gr".

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    estimate_g(sensor)
    for n, d in enumerate(["x", "y", "z"]):
        g_n = sensor.data["iss_g"][n]
        a = sensor.data["iss_a"+d]
        
        varname = "iss_a" + d + "_gr"
        sensor.data[varname] = a - g_n
        sensor.units[varname] = "m/s2"


def estimate_velocities(sensor):
    """ Use accelerations with g removed to estimate velocities in the iss.

    iss = initial sensor system
    Integrate the accelerations.

    Adds arrays containing the velocity: "iss_v{x,y,z}".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    dt = sensor.data["dt"]
    for d in ["x", "y", "z"]:
        a = sensor.data["iss_a" + d + "_gr"]
        v = np.cumsum(a*dt)

        varname = "iss_v" + d
        sensor.data[varname] = v
        sensor.units[varname] = "m/s"


def estimate_positions(sensor):
    """ Use velocities to estimate positions in the iss.

    iss = initial sensor system
    Integrate the velocities.

    Adds arrays containing the position: "iss_{x,y,z}".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    dt = sensor.data["dt"]
    for d in ["x", "y", "z"]:
        v = sensor.data["iss_v" + d]
        x = np.cumsum(v*dt)

        varname = "iss_" + d
        sensor.data[varname] = x
        sensor.units[varname] = "m"


def estimate_height(sensor):
    """ Estimate height in the lab frame.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    dt = sensor.data["dt"]
    a = sensor.data["lab_az_gr"]
    v = np.cumsum(a*dt)
    z = np.cumsum(v*dt)

    sensor.data["lab_vz"] = v
    sensor.units["lab_vz"] = "m/s"
    sensor.data["lab_z"] = z
    sensor.units["lab_z"] = "m"