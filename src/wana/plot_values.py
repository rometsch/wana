import os
import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import argparse

from pprint import pprint

from session import Session
from task import Task


def main():
    args = parse_cli_args()

    s = Session(args.session_file)

    pprint(s.start)
    pprint(s.stop)
    pprint((s.stop-s.start).seconds)

    t = Task(s, args.N_test)

    fl = t.sensors[0]
    fr = t.sensors[1]

    fl.integrate_angles()

    # pprint(t.meta)

    # exit(0)
    # fl.calc_g(2705300, 2705500)
    # fr.calc_g(2.731e6+160, 2.731e6+240)

    # fig = plot_acceleration_data_single([fl, fr])
    # fig = plot_gyro_data_single([fl, fr])
    fig = plot_angle_data_single(t.sensors)

    fig.suptitle(t.name)

    plt.show()


def plot_gyro_data(sensors):
    """ Plot gyroscope data.

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
    fig, axes = plt.subplots(3, N_sensors, sharex="all", sharey="row")
    sensor_axes = []
    if N_sensors == 1:
        sensor_axes = [axes]
    else:
        for n in range(N_sensors):
            sensor_axes.append(axes[:, n])

    #
    # Plot the data for all sensors
    #
    for axs, foot in zip(sensor_axes, sensors):
        for varname, ax in zip(["rx", "ry", "rz"], axs):
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")
            ax.grid(alpha=0.6)
            data_unit = foot.units[varname]
            ax.set_ylabel(f"{varname} [{data_unit}]")

        axs[0].set_title(foot.name)

        time_unit = foot.units["time"]
        axs[-1].set_xlabel(f"time [{time_unit}]")

    return fig


def plot_angle_data_single(sensors):
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


def plot_gyro_data_single(sensors):
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


def plot_acceleration_data_single(sensors):
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


def plot_acceleration_data(sensors):
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
    fig, axes = plt.subplots(4, N_sensors, sharex="all", sharey="row")
    sensor_axes = []
    if N_sensors == 1:
        sensor_axes = [axes]
    else:
        for n in range(N_sensors):
            sensor_axes.append(axes[:, n])

    #
    # Plot the data for all sensors
    #
    for axs, foot in zip(sensor_axes, sensors):
        for varname, ax in zip(["a", "ax", "ay", "az"], axs):
            x = foot.data["time"]
            y = foot.data[varname] / foot.g
            ax.plot(x, y, label=f"{varname}")
            ax.grid(alpha=0.6)
            data_unit = foot.units[varname]
            ax.set_ylabel(f"{varname} [{data_unit}]")

        axs[0].set_title(foot.name)

        time_unit = foot.units["time"]
        axs[-1].set_xlabel(f"time [{time_unit}]")

    return fig


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
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
