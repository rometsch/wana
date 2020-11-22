import argparse
import os

import numpy as np

from wana.session import Session
from wana.task import Task


def main():
    args = parse_cli_args()

    outdir = args.outdir

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.splitext(args.data_file)[-1] == ".zip":
        s = Sensor(args.data_file, "phyphox")
        sensors = [s]
        name = "phyphox"
        for s in sensors:
            write_data(s, s.name, outdir)

    else:
        s = Session(args.data_file)
        tasks = [Task(s, t) for t in range(len(s.tests))]
        for n, t in enumerate(tasks):
            for s in t.sensors:
                write_data(s, f"task_{n}_{s.name}", outdir)


def write_data(sen, name, outdir):
    """ Write data to a tab separated output file.

    Paramters
    ---------
    sen : 'obj`: sensor.Sensor
        Sensor to write out.
    name : str
        Name of the output file (excluding file extension).
    outdir : str
        Path to the output directory.
    """
    t = sen.data["time"]
    ax = sen.data["ax"]
    ay = sen.data["ay"]
    az = sen.data["az"]

    rx = sen.data["rx"]/180*np.pi
    ry = sen.data["ry"]/180*np.pi
    rz = sen.data["rz"]/180*np.pi

    header = "# Accelerometer and gyroscope data converted using wana."
    header += "\n# t = time in s"
    header += "\n# a{x,y,z} = acceleration in sensors coordinate system in m/s"
    header += "\n# r{x,y,z} = angular velocity in sensors coordinate system in rad/s"
    header += "\n# t\tax\tay\taz\trx\try\trz"

    filename = os.path.join(outdir, f"{name}.dat")

    output_format = "{:6e}"
    pattern = f"{output_format}\t" * 6 + output_format

    with open(filename, "w") as outfile:
        print(header, file=outfile)
        for n in range(len(t)):
            print(pattern.format(t[n], ax[n], ay[n], az[n], rx[n], ry[n], rz[n]),
                  file=outfile)


def parse_cli_args():
    """ Parse command line arguments.

    Returns
    -------
    Namespace
        A namespace containing the arguments as attributes.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_file", help="XML (gaitlab mobile) or zip (phyphox app) file.")
    parser.add_argument(
        "outdir", help="Output directory to write files to.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
