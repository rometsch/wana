import os

from wana.sensor import Sensor
from wana.calibration import load_calibration


class Task:
    """ Representation of a single test sequence. 

    Parameters
    ----------
    session: Session
        Session object containing the information.
    n: int
        Index of the test in the test list.
    """

    def __init__(self, session, n):
        self.session = session
        self.meta = session.tests[n]
        self.name = self.meta["Name"]
        self.desc = self.meta["ShortDescription"]
        self.data_dir = os.path.dirname(session.filename)
        self.sensors = []
        self.load_sensor_data()

    def load_sensor_data(self):
        for s in self.meta["MoteList"]["Mote"]:
            data_file = os.path.join(self.data_dir, s["File"])
            identifier = os.path.basename(data_file).split("_")[1]
            cal_acc, cal_gyro = load_calibration(identifier, self.data_dir)
            name = s["Position"]
            start_index = int(s["Tag"]["Start"])
            stop_index = int(s["Tag"]["Stop"])
            sensor = Sensor(
                data_file, f"{name}_{identifier}", cal_acc=cal_acc, cal_gyro=cal_gyro)
            sensor.set_time(0.01)
            sensor.trim_data(start_index, stop_index)
            sensor.integrate_angles()
            self.sensors.append(sensor)
