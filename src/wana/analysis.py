""" Functions to analyze the sensor data. """

import numpy as np
import numpy.ma as ma

# the numerical constant of earth's acceleration
g_constant = 9.81


def flag_resting(sensor, N_kernel=20):
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

    mask = delta < 0.02*g_constant

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
    print("{} steps detected for sensor {}".format(len(steps), sensor.name))


def regenerate_masks(sensor):
    """ Use the step intervals with padding to regenerate the resting and moving masks.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    N = len(sensor.data["counter"])
    mask_moving = np.zeros(N, dtype=bool)

    for I in sensor.data["interval_steps"]:
        mask_moving[I[0]:I[1]] = True

    mask_resting = np.logical_not(mask_moving)

    sensor.data["mask_resting"] = mask_resting
    sensor.data["mask_moving"] = mask_moving


def estimate_g(sensor, mode="average"):
    """ Estimate g vector from acceleration data in iss system.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    mode: str
        Select the algorithm to fit the data. Default: "average"
    """

    mode = "linear"
    extract_resting_accelerations(sensor)
    if mode == "average":
        estimate_g_average(sensor)
    if mode == "linear":
        estimate_g_poly(sensor, 1)
    set_g_units(sensor)


def extract_resting_accelerations(sensor):
    """ Extract accelerations for the times when sensor is at rest.

    iss = initial sensor system

    Adds masked arrays containing the accelerations when the sensor is at rest: "iss_a{x,y,z}_rest".


    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    mask_resting = np.logical_not(sensor.data["mask_resting"])
    sensor.data["time_resting"] = ma.masked_array(
        sensor.data["time"], mask=mask_resting)
    sensor.units["time_resting"] = "s"
    for d in ["x", "y", "z"]:
        a_masked = ma.masked_array(sensor.data["iss_a" + d], mask=mask_resting)

        varname = "iss_a" + d + "_rest"
        sensor.data[varname] = a_masked
        sensor.units[varname] = "m/s2"


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
    gx = np.average(sensor.data["iss_ax_rest"])
    gy = np.average(sensor.data["iss_ay_rest"])
    gz = np.average(sensor.data["iss_az_rest"])

    g = np.sqrt(gx**2 + gy**2 + gz**2)

    sensor.data["iss_g"] = np.array([gx, gy, gz])  # /g*g_constant
    sensor.units["iss_g"] = "m/s2"

    # add components
    N = len(sensor.data["counter"])
    sensor.data["iss_gx"] = np.ones(N)*gx
    sensor.data["iss_gy"] = np.ones(N)*gy
    sensor.data["iss_gz"] = np.ones(N)*gz


def estimate_g_poly(sensor, order):
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
    time_resting = sensor.data["time_resting"]
    time = sensor.data["time"]
    for d in ["x", "y", "z"]:
        acc = sensor.data[f"iss_a{d}_rest"]
        p = np.polyfit(time_resting, acc, order)
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
    mask_resting = sensor.data["mask_resting"]
    for n, d in enumerate(["x", "y", "z"]):
        varname = frame + "_a" + d + "_gr"
        sensor.data[varname][mask_resting] = 0

    calculate_norm(sensor, frame+"_a{}_gr", unit="m")


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
