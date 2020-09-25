import numpy as np
import numpy.ma as ma

# the numerical constant of earth's acceleration
g_constant = 9.81


def flag_resting(sensor, varname, N_kernel=20, order=1,
                 prior_mask=None, reference="max", delta_threshold=0.02):
    """ Get a mask for when the sensor is resting on ground.

    Compare the total acceleration to the value of g = 9.81 m/s2

    Adds a bool array "mask_resting.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    """
    values = sensor.data[varname]
    time = sensor.data["time"]
    if prior_mask is not None:
        fit_values = values[prior_mask]
        fit_time = time[prior_mask]
    else:
        fit_values = values
        fit_time = time

    if reference == "g":
        delta = np.abs(fit_values - g_constant)
        ref_val = g_constant
    elif reference == "max":
        p = np.polyfit(fit_time, fit_values, order)
        poly = np.poly1d(p)
        delta = np.abs(values - poly(time))
        ref_val = np.max(np.abs(values))
    else:
        raise ValueError("reference =", reference, "is not supported!")

    mask = delta < delta_threshold*ref_val

    if prior_mask is not None:
        mask = np.logical_and(mask, prior_mask)

    kernel = np.ones(N_kernel)/N_kernel
    conv = np.convolve(mask, kernel, mode="same")
    mask_resting = conv > 0.5
    mask_moving = conv <= 0.5

    sensor.data[f"mask_{varname}_resting"] = mask_resting
    sensor.data[f"mask_{varname}_moving"] = mask_moving

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
    resting = sensor.data["mask_a_resting"]

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


def find_step_intervals(sensor, min_time=0.5):
    """ Detect the bounding interval of steps.

    Adds 4 arrays:
        mask_start/stop_resting: boolean array with 1 when resting starts/stops
        index_start/stop_resting: integer array with indices when resting starts/stops

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    min_time: float
        Minimum time of a step.
    """
    time = sensor.data["time"]
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
        if time[up] - time[low] < min_time:
            continue
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

    sensor.data["mask_a_resting"] = mask_resting
    sensor.data["mask_a_moving"] = mask_moving
