""" Function to smooth the data. """
import numpy as np

def smooth_window(sensor, varname, N=5):
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
    
    kernel = np.ones(N)/N
    
    sensor.data[varname] = np.convolve(v,kernel)
    