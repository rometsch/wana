import numpy as np


def load_rawdata(datafile):
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
    return data_dict


class Sensor:
    """ Hold data about a foot. """

    def __init__(self, datafile, name, cal_acc=None, cal_gyro=None):
        self.datafile = datafile
        self.name = name
        self.data = load_rawdata(self.datafile)
        self.init_accelerations()
        self.init_unit_names()
        self.g = 1
        self.cal_gyro = cal_gyro
        self.cal_acc = cal_acc

        if cal_gyro is not None:
            self.calibrate_gyro()
        if cal_acc is not None:
            self.calibrate_acc()

    def init_accelerations(self):
        ax = self.data["ax"]
        ay = self.data["ay"]
        az = self.data["az"]
        self.data["a"] = np.sqrt(ax**2 + ay**2 + az**2)

    def init_unit_names(self):
        """ Initialize unit names for all data. """
        self.units = {"a": "raw",
                      "ax": "raw",
                      "ay": "raw",
                      "az": "raw",
                      "rx": "raw",
                      "ry": "raw",
                      "rz": "raw",
                      "angle_x": "raw",
                      "angle_y": "raw",
                      "angle_z": "raw",
                      "time": "s"}

    def set_time(self, timestep):
        """ Calculate a time variable using a start and stop value. 

        Parameters
        ----------
        timestep: float
            Time between two data points.
        """
        c = self.data["counter"]
        c = c - c[0]
        x = c/c[-1]
        self.data["time"] = c*timestep

    def trim_data(self, low, up):
        """ Trim the data to an interval defined by indices low and up.

        Parameters
        ----------
        low: int
            Start of the interval.
        up: int
            End of the interval.
        """
        print("total data length", len(self.data["ax"]))
        for key in self.data:
            self.data[key] = self.data[key][low:up+1]

    def calibrate_gyro(self):
        """ Use calibration for gyroscope values. """
        print(
            f"calibrating gyro for sensor {self.name} ({self.datafile}) using {self.cal_gyro.file_name}")
        offsets = self.cal_gyro.offset
        scaling = self.cal_gyro.scaling
        for varname, offset, scale in zip(["rx", "ry", "rz"], offsets, scaling):
            print(varname, "data = ", np.max(
                self.data[varname]), "offset = ", offset, "scaling = ", scale)
            self.data[varname] -= offset
            self.data[varname] /= scale
        self.units["rx"] = "deg/s"
        self.units["ry"] = "deg/s"
        self.units["rz"] = "deg/s"
        self.units["angle_x"] = "deg"
        self.units["angle_y"] = "deg"
        self.units["angle_z"] = "deg"

    def calibrate_acc(self):
        """ Use calibration for acceleration values. """
        print(
            f"calibrating acceleration for sensor {self.name} ({self.datafile}) using {self.cal_acc.file_name}")
        offsets = self.cal_acc.offset
        scaling = self.cal_acc.scaling
        for varname, offset, scale in zip(["ax", "ay", "az"], offsets, scaling):
            print(varname, "data = ", np.max(
                self.data[varname]), "offset = ", offset, "scaling = ", scale)
            self.data[varname] -= offset
            self.data[varname] /= scale
        self.units["ax"] = "9.81 m/s2"
        self.units["ay"] = "9.81 m/s2"
        self.units["ay"] = "9.81 m/s2"
        self.units["a"] = "9.81 m/s2"
        self.init_accelerations()

    def integrate_angles(self):
        """ Integrate angular velocity to get angles. """
        for d in ["x", "y", "z"]:
            # get angular_velocity in degree/seconds
            varname = "r" + d
            omega = self.data[varname]
            # get time in seconds
            time = self.data["time"]
            # dt = time[1] - time[0]
            dt = 0.01
            I = np.cumsum(omega)
            I *= dt
            self.data[f"angle_{d}"] = I
