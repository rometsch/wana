""" Functions to import data from the raw data files. """
import numpy as np


def load_rawdata_mobilegaitlab(datafile):
    """ Load data from the raw data file provided by the app.

    The data dict contains the following values:
        counter: time variable of the sensor
        a{x,y,z}: {x,y,z}-component of the accelerometer
        r{x,y,z}: {x,y,z}-component of the gyroscope
    The values are stored in the raw datafiles as
    int32 for the counter and int16 for the data values.


    Parameters
    ----------
    datafile: str
        Path to the datafile.

    Returns
    -------
    dict
        A dictionary containing the return values.
    """
    itype = np.int16
    rawdata = np.fromfile(datafile, dtype=[
        ("counter", np.int32),
        ("ax", itype),
        ("ay", itype),
        ("az", itype),
        ("rx", itype),
        ("ry", itype),
        ("rz", itype)])
    data_dict = {key: np.array(rawdata[key], dtype=float)
                 for key in rawdata.dtype.fields}
    N = len(data_dict["counter"])
    data_dict["dt"] = 0.01*np.ones(N, dtype=float)
    return data_dict
