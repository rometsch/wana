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

    if args.list_vars:
        print("Available sensor data:")
        for key in t.sensors[0].data:
            print(key)

    for key in args.plots:
        if key == "manual":
            fig = plot_functions[key](
                t.sensors, args.names, show_steps=args.show_steps)
        else:
            fig = plot_functions[key](t.sensors, show_steps=args.show_steps)
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
    parser.add_argument("--show-steps", default=False, action="store_true",
                        help="Show the steps as shaded intervals.")
    parser.add_argument("--list-vars", default=False, action="store_true",
                        help="Print a list of available variables.")
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


def plot_vars(sensors, varnames, y_name="", show_steps=False):
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
        units = []
        for varname in varnames:
            x = foot.data["time"]
            y = foot.data[varname]
            ax.plot(x, y, label=f"{varname}")
            try:
                units.append(foot.units[varname])
            except KeyError:
                units.append(None)

        ax.set_title(foot.name)

        ax.grid(alpha=0.6)
        ax.legend()

        time_unit = foot.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        if all([units[0] == u for u in units]) and units[0] is not None:
            data_unit = units[0]
            ax.set_ylabel(f"{y_name} [{data_unit}]")
        else:
            ax.set_ylabel(f"arbitrary units")

        try:
            if show_steps:
                steps = foot.data["interval_steps"]
                for step in steps:
                    t = foot.data["time"]
                    t_shade = t[step[0]:step[1]]
                    ax.fill_between(t_shade, 1, alpha=0.1, color="k",
                                    transform=ax.get_xaxis_transform())
        except KeyError:
            pass

    return fig


def plot_angles(sensors, **kwargs):
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
    return plot_vars(sensors, ["angle_x", "angle_y", "angle_z"], y_name="angle", **kwargs)


def plot_gyro(sensors, **kwargs):
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
    return plot_vars(sensors, ["rx", "ry", "rz"], y_name="$\omega$", **kwargs)


def plot_acceleration(sensors, **kwargs):
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
    return plot_vars(sensors, ["a", "ax", "ay", "az"], y_name="a", **kwargs)


def plot_lab_3d_trajectory(sensors, show_steps=False):
    """ Plot the trajectory in 3d in the lab frame.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.

    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    for sensor in sensors:
        x = sensor.data["lab_x_lc"]
        y = sensor.data["lab_y_lc"]
        z = sensor.data["lab_z_lc"]
        line, = ax.plot(x, y, z, label=sensor.name)
        color = line.get_color()
        if show_steps:
            for I in sensor.data["interval_steps"]:
                x_i = x[I[0]]
                y_i = y[I[0]]
                z_i = z[I[0]]
                ax.scatter(x_i, y_i, z_i, marker="d", color=color)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    axisEqual3D(ax)

    ax.legend()

    return fig


def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))()
                        for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


plot_functions = {"acc": plot_acceleration,
                  "angle": plot_angles,
                  "rot": plot_gyro,
                  "manual": plot_vars,
                  "trajectory": plot_lab_3d_trajectory
                  }
available_plots = [key for key in plot_functions]


if __name__ == "__main__":
    main()
