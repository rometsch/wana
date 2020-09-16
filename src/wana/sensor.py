import numpy as np
import wana.analysis as analysis
import wana.transformations as trafo


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
    N = len(data_dict["counter"])
    data_dict["dt"] = 0.01*np.ones(N, dtype=float)
    return data_dict


class Sensor:
    """ Hold data about a foot. """

    def __init__(self, datafile, name, sample_rate=None, cal_acc=None, cal_gyro=None, trim_low=None, trim_up=None):
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
        if sample_rate is not None:
            self.set_time(1/sample_rate)
            
        if trim_low is not None and trim_up is not None:
            self.trim_data(trim_low, trim_up)

        if all([x is not None for x in [sample_rate, cal_acc, cal_gyro]]):
            self.postprocess()

    def postprocess(self):
        self.integrate_angles()
        self.calc_delta_angle()
        trafo.transform_to_reference_system(self)
        analysis.flag_resting(self)
        analysis.estimate_g(self)
        trafo.calc_lab_ez(self)
        trafo.calc_lab_ehor(self)
        trafo.calc_rotation_iss_to_lab(self)
        analysis.remove_g(self)
        analysis.estimate_velocities(self, "iss")
        analysis.estimate_positions(self, "iss")
        trafo.iss_to_lab(self)
        analysis.estimate_velocities(self, "lab")
        analysis.estimate_positions(self, "lab")
        try:
            analysis.find_step_intervals(self)
            analysis.estimate_velocities(self, "iss", perstep=True)
            analysis.estimate_positions(self, "iss", perstep=True)

            analysis.estimate_velocities(self, "lab", perstep=True)
            analysis.estimate_positions(self, "lab", perstep=True)
        except IndexError:
            print("Could not detect any steps!")

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
            self.data[varname] *= 9.81/scale
        self.units["ax"] = "m/s2"
        self.units["ay"] = "m/s2"
        self.units["az"] = "m/s2"
        self.units["a"] = "m/s2"
        self.init_accelerations()

    def calc_delta_angle(self):
        """Calculate the change of angle for the timesteps.

        Parameters
        ----------
        dt: float
            Time between datapoints in seconds.
        """
        for d in ["x", "y", "z"]:
            varname = "delta_angle_" + d
            omega = self.data["r"+d]
            dt = self.data["dt"]
            delta = omega*dt
            self.data[varname] = delta
            self.units[varname] = "rad"

    def integrate_angles(self):
        """ Integrate angular velocity to get angles. 

        Parameters
        ----------
        dt: float
            Time between datapoints in seconds.
        """
        for d in ["x", "y", "z"]:
            # get angular_velocity in degree/seconds
            varname = "r" + d
            omega = self.data[varname]
            # get time in seconds
            dt = self.data["dt"]
            I = np.cumsum(omega)
            I *= dt

            varname = f"angle_{d}"
            self.data[varname] = I
            self.units[varname] = "rad"
