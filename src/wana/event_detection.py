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
