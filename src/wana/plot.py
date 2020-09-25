import argparse
import sys
import os

import matplotlib.pyplot as plt
import numpy as np

from wana.sensor import Sensor
from wana.session import Session
from wana.task import Task


def main():
    args = parse_cli_args()

    if os.path.splitext(args.data_file)[-1] == ".zip":
        s = Sensor(args.data_file, "phyphox")
        sensors = [s]
        name = "phyphox"
    else:
        s = Session(args.data_file)
        t = Task(s, args.task)
        sensors = t.sensors
        name = t.name
        if args.sensor is not None:
            sensors = [sensors[args.sensor]]

    if args.smooth is not None:
        for s in sensors:
            s.smooth(args.smooth)
            s.postprocess()
    if args.trim is not None:
        for s in sensors:
            low = args.trim[0]
            up = args.trim[1]
            s.trim_time(low, up)

    if args.list_vars:
        print("Available sensor data:")
        for key in sensors[0].data:
            print(key)

    # plot_acceleration([s, s])

    # plt.show()

    for key in args.plots:
        if key == "manual":
            fig = plot_functions[key](
                sensors, args.names, show_steps=args.show_steps, stepwise=args.stepwise)
        elif key in available_plots:
            fig = plot_functions[key](
                sensors, show_steps=args.show_steps, stepwise=args.stepwise)
        else:
            fig = plot_vars(
                sensors, [key], show_steps=args.show_steps, stepwise=args.stepwise)

        fig.suptitle(name)

        add_footnote(fig, "wana plot " + " ".join(sys.argv[1:]))
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
        "data_file", help="XML (gaitlab mobile) or zip (phyphox app) file.")
    parser.add_argument("-p", "--plots", nargs="+", type=str,
                        default=[k for k in available_plots if k not in [
                            "manual", "trajectory"]],
                        help="Plots to produce.")
    parser.add_argument("-o", "--outfile", help="File to output data to.")
    parser.add_argument("-n", "--names", nargs="+",
                        type=str, help="Variable names to plot.")
    parser.add_argument("--show-steps", default=False, action="store_true",
                        help="Show the steps as shaded intervals.")
    parser.add_argument("--list-vars", default=False, action="store_true",
                        help="Print a list of available variables.")
    parser.add_argument("-s", "--stepwise", default=False, action="store_true",
                        help="Use data based on single steps.")
    parser.add_argument("--trim", nargs=2, type=float,
                        help="Trim data to time.")
    parser.add_argument("--smooth", type=int,
                        help="Smooth over N datapoints.")
    parser.add_argument("--task", type=int, default=0,
                        help="For gaitlab files: select the task.")
    parser.add_argument("--sensor", type=int,
                        help="Select a single sensor to display.")
    args = parser.parse_args()

    if "manual" in args.plots and args.names is None:
        print("For the manual plot, at least one variables name must be specified with -n/--names!")
        parser.print_help()
        exit(1)

    return args


def add_footnote(fig, text):
    """ Add a footnote to a plot.

    Parameters
    ----------
    fig: pyplot.figure
        Figure to put the footnote to.
    text: str
        String to put into the footnote.
    """
    fig.text(0.02, 0.02, text, fontdict={"fontsize": 4})


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


def expand_varnames(varnames, vector_pattern="VEC"):
    """ Expand variable names to vectors if indicated.

    Vector quantities are indicated by placing the vector_pattern string
    at the position where x,y,z would be placed.

    I.e. with vector_pattern="VEC"
    iss_aVEC_gr -> iss_ax_gr, iss_ay_gr, iss_az_gr

    Parameters
    ----------
    varnames: list of str
        Variable names to be expanded.
    vector_pattern: str
        Pattern to be replaced.

    Returns
    -------
    list of str
        Expanded variables names.
    """
    to_expand = [n for n in varnames if vector_pattern in n]
    to_leave = [n for n in varnames if not n in to_expand]
    expanded = []
    for n in to_expand:
        for d in ["x", "y", "z"]:
            new = n.replace(vector_pattern, d)
            expanded.append(new)
    rv = to_leave + expanded
    return rv


def plot_vars(sensors, varnames, y_name="", show_steps=False, stepwise=False):
    """ Plot angle data.

    Parameters
    ----------
    sensors: list of Sensor
        A list of the sensors to be plotted.
    varnames: list of str
        Variable names to be plotted.
    show_steps: bool
        Indicated detected steps as shaded region.
    stepwise: bool
        Use data based on single steps rather than full length.
    Returns
    -------
    plt.fig
        Pyplot figure holding the plot.
    """
    varnames = expand_varnames(varnames)

    N_sensors = len(sensors)
    fig, axes = plt.subplots(1, N_sensors, sharex="all", sharey="row")
    if N_sensors == 1:
        axes = [axes]

    if stepwise:
        steps = []
        steps_axes = []
        for sensor, ax in zip(sensors, axes):
            for step in sensor.steps:
                steps.append(step)
                steps_axes.append(ax)
        plot_sensors = steps
        plot_axes = steps_axes
    else:
        plot_axes = axes
        plot_sensors = sensors

    colors = [c["color"] for c in plt.rcParams['axes.prop_cycle']]

    for sensor, ax in zip(plot_sensors, plot_axes):

        units = []
        for varname, color in zip(varnames, colors):
            x = sensor.data["time"]
            y = sensor.data[varname]
            ax.plot(x, y, color=color, alpha=1)
            try:
                units.append(sensor.units[varname])
            except KeyError:
                units.append(None)

        ax.set_title(sensor.name)

        ax.grid(alpha=0.6)

        time_unit = sensor.units["time"]
        ax.set_xlabel(f"time [{time_unit}]")

        if all([units[0] == u for u in units]) and units[0] is not None:
            data_unit = units[0]
            ax.set_ylabel(f"{y_name} [{data_unit}]")
        else:
            ax.set_ylabel(f"arbitrary units")

        try:
            if show_steps:
                steps = sensor.data["interval_steps"]
                for step in steps:
                    t = sensor.data["time"]
                    t_shade = t[step[0]:step[1]]
                    ax.fill_between(t_shade, 1, alpha=0.1, color="k",
                                    transform=ax.get_xaxis_transform())
        except KeyError:
            pass

    for ax in axes:
        labels = varnames
        lines = ax.lines[:len(varnames)]
        ax.legend(lines, labels)

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


def plot_lab_3d_trajectory(sensors, show_steps=False, **kwargs):
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
