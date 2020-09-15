import argparse
import os
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

from wana.session import Session
from wana.task import Task


def main():
    args = parse_cli_args()

    s = Session(args.session_file)
    t = Task(s, args.N_test)

    for key in args.plots:
        if key == "manual":
            fig = plot_functions[key](t.sensors, args.names)
        else:
            fig = plot_functions[key](t.sensors)
        fig.suptitle(t.name)
        if args.outfile:
            append_name = len(args.plots) > 1
            save_plot(fig, key, args.outfile, append_name=append_name)

    if args.outfile is None:
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
    parser.add_argument("-o", "--outfile", help="File to output data to.")
    parser.add_argument("-n", "--names", nargs="+",
                        type=str, help="Variable names to plot.")
    args = parser.parse_args()

    if "manual" in args.plots and args.names is None:
        print("For the manual plot, at least one variables name must be specified with -n/--names!")
        parser.print_help()
        exit(1)

    return args


def save_plot(fig, name, filename, append_name=False):
    """ Save a plot as image file.

    The type is determined by the filename ending.

    Paramters
    ---------
    fig: pyplot.figure
        The figure object holding the plot.
    name: str
        Name of the plot.
    filename: str
        Path of the output image.
    append_name: bool
        Whether or not to append the plot name before the extension.
    """
    if append_name:
        filename = append_before_extension(filename, name)
    fig.savefig(filename, dpi=300)


def append_before_extension(s, i):
    """ Insert i before the extension of s. 

    Example: s = foo.bar, i = baz, result = foo_baz.bar

    Parameters
    ----------
    s: str
        Base string.
    i: str
        Part to be inserted before last '.'.         
    """
    parts = s.split(".")
    if len(parts) == 1:
        rv = s + "_" + i
    else:
        rv = ".".join(parts[:-1]) + "_" + i + "." + parts[-1]
    return rv


def plot_vars(sensors, varnames):
    """ Plot angle data.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.
    varnames: list of str
        Variable names to be plotted.

    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """
    N_sensors = len(sensors)
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")

    for foot, ax in zip(sensors, axes):
        for varname in varnames:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")

        ax.set_title(foot.name)

        ax.grid(alpha=0.6)
        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        ax.set_ylabel(f"in arbitrary units")

    return fig


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
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")

    for foot, ax in zip(sensors, axes):
        for varname in ["angle_x", "angle_y", "angle_z"]:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")

        ax.set_title(foot.name)

        ax.grid(alpha=0.6)
        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        data_unit = foot.units["angle_x"]
        ax.set_ylabel(f"angle [{data_unit}]")

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

    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")

    for foot, ax in zip(sensors, axes):
        for varname in ["rx", "ry", "rz"]:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")

        ax.set_title(foot.name)

        ax.grid(alpha=0.6)
        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        data_unit = foot.units["rx"]
        ax.set_ylabel(f"$\omega$ [{data_unit}]")

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

    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")

    for foot, ax in zip(sensors, axes):
        for varname in ["a", "ax", "ay", "az"]:
            x = foot.data["time"]
            y = foot.data[varname] / foot.g
            ax.plot(x, y, label=f"{varname}")

        ax.set_title(foot.name)

        ax.grid(alpha=0.6)
        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        data_unit = foot.units["a"]
        ax.set_ylabel(f"a [{data_unit}]")

    return fig


plot_functions = {"acc": plot_acceleration,
                  "angle": plot_angles,
                  "rot": plot_gyro,
                  "manual": plot_vars}
available_plots = [key for key in plot_functions]


if __name__ == "__main__":
    main()
