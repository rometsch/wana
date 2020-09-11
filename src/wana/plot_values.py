import argparse
import os
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

from session import Session
from task import Task


def main():
    args = parse_cli_args()

    s = Session(args.session_file)
    t = Task(s, args.N_test)

    for key in args.plots:
        fig = plot_functions[key](t.sensors)
        fig.suptitle(t.name)

    plt.show()


def parse_cli_args():
    """ Parse command line arguments. 

    Returns
    -------
    Namespace
        A namespace containing the arguments as attributes.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "session_file", help="XML file describing the session.")
    parser.add_argument("N_test", type=int, help="Number of test to load.")
    parser.add_argument("-p", "--plots", nargs="+", type=str, choices=available_plots,
                        default=available_plots, help="Plots to produce.")
    args = parser.parse_args()
    return args


def plot_angles(sensors):
    """ Plot angle data.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.

    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """
    N_sensors = len(sensors)
    #
    # Make a grid of coordinate systems
    #
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")
    #
    # Plot the data for all sensors
    #
    for foot, ax in zip(sensors, axes):
        for varname in ["angle_x", "angle_y", "angle_z"]:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")
            ax.grid(alpha=0.6)
            data_unit = foot.units[varname]
            ax.set_ylabel(f"{varname} [{data_unit}]")

        ax.set_title(foot.name)

        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

    return fig


def plot_gyro(sensors):
    """ Plot gyro data.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.

    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """
    N_sensors = len(sensors)
    #
    # Make a grid of coordinate systems
    #
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")
    #
    # Plot the data for all sensors
    #
    for foot, ax in zip(sensors, axes):
        for varname in ["rx", "ry", "rz"]:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")
            ax.grid(alpha=0.6)
            data_unit = foot.units[varname]
            ax.set_ylabel(f"{varname} [{data_unit}]")

        ax.set_title(foot.name)

        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

    return fig


def plot_acceleration(sensors):
    """ Plot acceleration data.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.

    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """
    N_sensors = len(sensors)
    #
    # Make a grid of coordinate systems
    #
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")
    #
    # Plot the data for all sensors
    #
    for foot, ax in zip(sensors, axes):
        for varname in ["a", "ax", "ay", "az"]:
            x = foot.data["time"]
            y = foot.data[varname] / foot.g
            ax.plot(x, y, label=f"{varname}")
            ax.grid(alpha=0.6)
            data_unit = foot.units[varname]
            ax.set_ylabel(f"{varname} [{data_unit}]")

        ax.set_title(foot.name)

        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

    return fig


plot_functions = {"acc": plot_acceleration,
                  "angle": plot_angles,
                  "rot": plot_gyro}
available_plots = [key for key in plot_functions]


if __name__ == "__main__":
    main()
