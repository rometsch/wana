import os

import numpy as np
import wana.analysis as analysis
import wana.transformations as trafo
import wana.smooth as smooth
import wana.load as load
import wana.event_detection as ed
import wana.integration as integration


class Sensor:
    """ Hold data about a foot. """

    def __init__(self, datafile, name, sample_rate=None, cal_acc=None, cal_gyro=None, trim_low=None, trim_up=None):
        self.datafile = datafile
        self.name = name
        if os.path.splitext(datafile)[-1] == ".zip":
            self.data = load.load_rawdata_phyphox(self.datafile)
        else:
            self.data = load.load_rawdata_mobilegaitlab(self.datafile)
        self.init_accelerations()
        self.init_unit_names()
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

        # if all([x is not None for x in [sample_rate, cal_acc, cal_gyro]]):
        # try:
        #     self.postprocess()
        # except TypeError:
        #     pass
        self.postprocess()

        self.steps = []


    def postprocess(self):
        self.integrate_angles()
        self.calc_delta_angle()

        trafo.transform_to_reference_system(self)

        ed.flag_resting(self, "a", reference="g", delta_threshold=0.02)
        analysis.extract_resting_values(self, "a")

        analysis.estimate_g(self, order=1)

        trafo.calc_lab_ez(self)
        trafo.calc_lab_ehor(self)

        trafo.calc_trafo_iss_to_lab(self)
        trafo.iss_to_lab(self, "a{}", unit="m/s2")
        # trafo.iss_to_lab_acc(self)

        analysis.remove_g(self)
        trafo.iss_to_lab(self, "a{}_gr", unit="m/s2")

        integration.estimate_velocities(self, "iss")
        integration.estimate_positions(self, "iss", "v{}")

        integration.estimate_velocities(self, "lab")
        integration.estimate_positions(self, "lab", "v{}")
        try:
            ed.find_step_intervals(self)

            #
            # iss frame
            #
            integration.estimate_velocities(self, "iss", perstep=True)

            for name in ["iss_vx", "iss_vy", "iss_vz"]:
                integration.linear_step_correction(self, name, unit="m/s")

            integration.estimate_positions(self, "iss", "v{}_lc")
            integration.estimate_positions(self, "iss", "v{}_step")

            #
            # lab frame
            #
            integration.estimate_velocities(self, "lab", "v{}")
            for name in ["lab_vx", "lab_vy", "lab_vz"]:
                integration.linear_step_correction(self, name, unit="m/s")

            integration.estimate_positions(self, "lab", "v{}_lc")
            integration.estimate_positions(self, "lab", "v{}_lc", perstep=True)

            integration.estimate_positions(self, "lab", "v{}_step")
            analysis.calculate_velocity_direction_angles(self, "lab_v{}_lc")

        except IndexError:
            print("Could not detect any steps!")

    def smooth(self, N):
        """ Smooth the sensor data. """
        for var in ["ax", "ay", "az", "rx", "ry", "rz"]:
            smooth.smooth_window(self, var, N=N)

    def init_accelerations(self):
        ax = self.data["ax"]
        ay = self.data["ay"]
        az = self.data["az"]
        self.data["a"] = np.sqrt(ax**2 + ay**2 + az**2)

    def init_unit_names(self):
        """ Initialize unit names for all data. """
        self.units = {"a": "m/s",
                      "ax": "m/s",
                      "ay": "m/s",
                      "az": "m/s",
                      "rx": "deg/s",
                      "ry": "deg/s",
                      "rz": "deg/s",
                      "angle_x": "deg",
                      "angle_y": "deg",
                      "angle_z": "deg",
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

    def trim_time(self, low, up):
        """ Trim the data to an interval defined by indices low and up.

        Parameters
        ----------
        low: int
            Start of the interval.
        up: int
            End of the interval.
        """
        print("total data length", len(self.data["ax"]))
        time = self.data["time"]
        mask = np.logical_and(time >= low, time <= up)
        for key in self.data:
            if not isinstance(self.data[key], np.ndarray) or not len(self.data[key]) == len(mask):
                continue
            self.data[key] = self.data[key][mask]
        self.postprocess()

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

    def generate_steps(self):
        """ Construct an array with sensor objects for each step. """
        for n in range(1, len(self.data["interval_steps"]) + 1):
            self.generate_step(n)

    def generate_step(self, n):
        """ Construct a single step. """
        intervals = [[0, 0]] + self.data["interval_steps"] + [[-1, -1]]
        try:
            low = intervals[n-1][1]
            up = intervals[n+1][0]
            s = Step(self, low, up,
                     name=self.name + f"_step_{n}")
            self.steps.append(s)
        except (TypeError, IndexError):
            print(
                f"Something went wrong in analysing step {n} of sensor {self.name}!")


class Step(Sensor):

    def __init__(self, other, low, up, name=None):
        self.datafile = other.datafile
        if name is None:
            self.name = other.name
        else:
            self.name = name
        self.cal_gyro = other.cal_gyro
        self.cal_acc = other.cal_acc
        self.units = other.units

        self.data = {}
        for var in ["ax", "ay", "az", "rx", "ry", "rz", "counter", "dt", "time"]:
            self.data[var] = other.data[var][low:up]

        self.init_accelerations()

        self.postprocess()
