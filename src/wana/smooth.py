""" Function to smooth the data. """
import numpy as np
from scipy import signal


def smooth_window(sensor, varname, N=10):
    """ Run a moving average over the data.

    Parameres
    ---------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    varname: str
        Variable to smooth.
    N: int
        Length of the smoothing interval.
    """
    v = sensor.data[varname]

    N = int(N)
    kernel = np.ones(N)/N

    sensor.data[varname] = np.convolve(v, kernel, mode="same")


def smooth_lowpass(sensor, varname, fcrit=10):
    """ Run a lowpass filter over the data.

    Parameres
    ---------
    sensor: wana.sensor.Sensor
        Sensor object holding the data.
    varname: str
        Variable to smooth.
    """
    dt = sensor.data["dt"][0]
    sig = sensor.data[varname]
    sampling_rate = 1/dt
    filter_order = 10
    sos = signal.butter(filter_order, fcrit, 'lp', fs=sampling_rate, output='sos')
    filtered = signal.sosfilt(sos, sig)
    sensor.data[varname] = filtered
