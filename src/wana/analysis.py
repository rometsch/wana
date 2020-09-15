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

    N_kernel = 20
    kernel = np.ones(N_kernel)/N_kernel
    conv = np.convolve(mask, kernel, mode="same")
    sensor.data["mask_resting"] = conv > 0.5
    sensor.data["mask_moving"] = conv <= 0.5

    find_resting_edge(sensor)


def find_resting_edge(sensor):
    """ Detect the beginning of a resting interval.

    Adds 4 arrays:
        mask_start/stop_resting: boolean array with 1 when resting starts/stops
        index_start/stop_resting: integer array with indices when resting starts/stops

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    resting = sensor.data["mask_resting"]

    trigger_val = 0.5
    mask1 = (resting[:-1] < trigger_val) & (resting[1:] > trigger_val)
    mask2 = (resting[:-1] > trigger_val) & (resting[1:] < trigger_val)

    index_start_resting = np.flatnonzero(mask1)+1
    mask_start_resting = np.zeros(len(resting))
    mask_start_resting[index_start_resting] = 1

    sensor.data["index_start_resting"] = index_start_resting
    sensor.data["mask_start_resting"] = mask_start_resting

    index_stop_resting = np.flatnonzero(mask2)+1
    mask_stop_resting = np.zeros(len(resting))
    mask_stop_resting[index_stop_resting] = 1

    sensor.data["index_stop_resting"] = index_stop_resting
    sensor.data["mask_stop_resting"] = mask_stop_resting


def find_step_intervals(sensor):
    """ Detect the bounding interval of steps.

    Adds 4 arrays:
        mask_start/stop_resting: boolean array with 1 when resting starts/stops
        index_start/stop_resting: integer array with indices when resting starts/stops

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    steps = []
    index_stop = sensor.data["index_start_resting"]
    index_start = sensor.data["index_stop_resting"]

    if index_start[0] > index_stop[0]:
        index_stop = index_stop[1:]

    N_start = len(index_start)
    N_stop = len(index_stop)

    N_steps = min(N_start, N_stop)

    for n in range(N_steps):
        low = index_start[n]
        up = index_stop[n]
        l = (up-low)
        pad = int(0.1*l)
        steps.append([low-pad, up+pad])


    sensor.data["interval_steps"] = np.array(steps)


def estimate_g(sensor):
    """ Estimate g vector from acceleration data in iss system.

    iss = initial sensor system

    Adds masked arrays containing the accelerations when the sensor is at rest: "iss_a{x,y,z}_rest".

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
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

    sensor.data["iss_g"] = np.array([gx, gy, gz])#/g*g_constant
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
    a = sensor.data["lab_az"]
    v = np.cumsum(a*dt)
    z = np.cumsum(v*dt)

    sensor.data["lab_vz"] = v
    sensor.units["lab_vz"] = "m/s"
    sensor.data["lab_z"] = z
    sensor.units["lab_z"] = "m"


def estimate_height_step(sensor):
    """ Estimate height in the lab frame individually per step.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    dt = sensor.data["dt"]
    a = sensor.data["lab_az"]

    index_step_start = sensor.data["interval_steps"][:,0]
    # index_step_start = sensor.data["index_stop_resting"]

    v = np.cumsum(a*dt)
    for i in index_step_start:
        v[i:] -= v[i]

    z = np.cumsum(v*dt)

    for i in index_step_start:
        z[i:] -= z[i]


    sensor.data["lab_vz_step"] = v
    sensor.units["lab_vz_step"] = "m/s"
    sensor.data["lab_z_step"] = z
    sensor.units["lab_z_step"] = "m"
