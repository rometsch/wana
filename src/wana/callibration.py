import numpy as np



def load_callibration(name, basedir="."):
    cal_acc = Callibration(basedir + "/CalibrationFiles/"
                           + name + "/2019-12-12_09-30/"
                           + name + "_acc.csv")
    cal_gyro = Callibration(basedir + "/CalibrationFiles/"
                            + name + "/2019-12-12_09-30/"
                            + name + "_gyro.csv")
    return cal_acc, cal_gyro


class Callibration:
    """ Load a GaitLab callibration file. """

    def __init__(self, file_name):
        self.file_name = file_name
        self.load()

    def load(self):
        self.data = np.genfromtxt(self.file_name, delimiter=",")
        self.offset = self.data[:, 0]
        self.crosstalk = self.data[:, 1:4]
        self.scaling_matrix = self.data[:, 4:7]
        sm = self.scaling_matrix
        self.scaling = [sm[0, 0], sm[1, 1], sm[2, 2]]
