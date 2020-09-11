import os
import sys

import numpy as np


def load_callibration(name, data_dir):
    try:
        cal_acc = Callibration(os.path.join(data_dir, name + "_acc.csv"))
    except OSError:
        print(f"Could not find acceleration callibration file for sensor {name}!"
              + "\nPlease copy it to the data directory.",
              file=sys.stderr)
        cal_acc = None
    try:
        cal_gyro = Callibration(os.path.join(data_dir, name + "_gyro.csv"))
    except OSError:
        print(f"Could not find gyroscope callibration file for sensor {name}!"
              + "\nPlease copy it to the data directory.",
              file=sys.stderr)
        cal_gyro = None
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
