import os
import sys

import numpy as np


def load_calibration(name, data_dir):
    try:
        cal_acc = Calibration(os.path.join(data_dir, name + "_acc.csv"))
    except OSError:
        print(f"Could not find acceleration calibration file for sensor {name}!"
              + "\nPlease copy it to the data directory.",
              file=sys.stderr)
        cal_acc = None
    try:
        cal_gyro = Calibration(os.path.join(data_dir, name + "_gyro.csv"))
    except OSError:
        print(f"Could not find gyroscope calibration file for sensor {name}!"
              + "\nPlease copy it to the data directory.",
              file=sys.stderr)
        cal_gyro = None
    return cal_acc, cal_gyro


class Calibration:
    """ Load a GaitLab calibration file. """

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

    def apply(self, vec):
        offsets = self.offset
        scaling = self.scaling
        # crosstalk_correction = np.linalg.inv(self.crosstalk)
        crosstalk_correction = self.crosstalk
        # vec = np.dot(crosstalk_correction, vec)
        for n, offset, scale in zip(range(3), offsets, scaling):
            vec[n] += -offset
            vec[n] /= scale
        return vec
