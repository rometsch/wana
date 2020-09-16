import os

from wana.calibration import load_calibration
from wana.sensor import Sensor


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
            sampling_rate = float(s["SamplingRate"])
            sensor = Sensor(
                data_file, f"{name}_{identifier}", 
                cal_acc=cal_acc, 
                cal_gyro=cal_gyro, 
                sample_rate=sampling_rate, 
                trim_low=start_index, 
                trim_up=stop_index)

            self.sensors.append(sensor)
