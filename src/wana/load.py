""" Functions to import data from the raw data files. """
import numpy as np
import zipfile

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


def load_rawdata_phyphox(datafile):
    """ Load data from the raw data file provided by the phyphox app.

    The data dict contains the following values:
        counter: time variable of the sensor
        a{x,y,z}: {x,y,z}-component of the accelerometer

    Parameters
    ----------
    datafile: str
        Path to the datafile.

    Returns
    -------
    dict
        A dictionary containing the return values.
    """
    with zipfile.ZipFile(datafile) as datazip:
        with datazip.open('Accelerometer.csv', "r") as infile:
            rawdata = infile.readlines()
            data = {}
            data["time"] = np.genfromtxt(rawdata, usecols=0, skip_header=1)
            N_acc = len(data["time"])
            data["ax"] = np.genfromtxt(rawdata, usecols=1, skip_header=1)
            data["ay"] = np.genfromtxt(rawdata, usecols=2, skip_header=1)
            data["az"] = np.genfromtxt(rawdata, usecols=3, skip_header=1)
            data["counter"] = np.arange(N_acc)
            data["dt"] = np.ones(N_acc) * (data["time"][1] - data["time"][0])
        with datazip.open('Gyroscope.csv', "r") as infile:
            rawdata = infile.readlines()
            data["rx"] = np.genfromtxt(rawdata, usecols=1, skip_header=1)
            data["ry"] = np.genfromtxt(rawdata, usecols=2, skip_header=1)
            data["rz"] = np.genfromtxt(rawdata, usecols=3, skip_header=1)
            N_gyro = len(data["rx"])
        N = min(N_acc, N_gyro)
        for key in data:
            data[key] = data[key][:N]
    return data
