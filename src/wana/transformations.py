""" Transformations of coordinate systems. """
import numpy as np
from scipy.spatial.transform import Rotation as R
from wana import analysis


def transform_to_reference_system(sensor):
    """ Tranform the data in the given sensor object to a reference system. 

    Use the angular velocities and time intervals to calculate the angles
    and then rotate at each time step to get the accelerations in the reference system.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Sensor to be transformed.
    """
    construct_rotations(sensor)
    transform_accelerations(sensor)


def construct_rotations(sensor):
    """ Construct rotations for every time step to transform to the inertial system.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the gyroscope data.
    """

    N = len(sensor.data["delta_angle_x"])

    rotations = [R.identity()]

    for n in range(1, N):
        # get angles in radians from the previous interval
        phi_x = sensor.data["delta_angle_x"][n-1]*np.pi/180
        phi_y = sensor.data["delta_angle_y"][n-1]*np.pi/180
        phi_z = sensor.data["delta_angle_z"][n-1]*np.pi/180

        # calculate the SORA rotation vector
        v = np.array([phi_x, phi_y, phi_z])
        r_step = R.from_rotvec(-v)

        r_prev = rotations[-1]

        # print("r_step", r_step.as_matrix())
        # print("r_prev", r_prev.as_matrix())

        r = r_step*r_prev

        rotations.append(r)

    sensor.data["rotation_to_iss"] = rotations


def transform_accelerations(sensor):
    """ Rotate acceleration vectors to the iss system.

    iss = initial sensor system

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the sensor data.
    """
    N = len(sensor.data["delta_angle_x"])

    sensor.data["iss_ax"] = np.ones(N)
    sensor.data["iss_ay"] = np.ones(N)
    sensor.data["iss_az"] = np.ones(N)

    sensor.units["iss_ax"] = "m/s2"
    sensor.units["iss_ay"] = "m/s2"
    sensor.units["iss_az"] = "m/s2"

    for n in range(0, N):
        r = sensor.data["rotation_to_iss"][n]
        vec = np.array([
            sensor.data["ax"][n],
            sensor.data["ay"][n],
            sensor.data["az"][n]
        ])

        vec_rot = r.apply(vec)

        sensor.data["iss_ax"][n] = vec_rot[0]
        sensor.data["iss_ay"][n] = vec_rot[1]
        sensor.data["iss_az"][n] = vec_rot[2]


def calc_lab_ez(sensor):
    """ Calculate the vertical unit vector of the lab frame.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the sensor data.
    """

    g_vec = np.array([
        sensor.data["iss_gx"][0],
        sensor.data["iss_gy"][0],
        sensor.data["iss_gz"][0]
        ])

    e_z = g_vec / np.linalg.norm(g_vec)

    varname = "lab_ez"
    sensor.data[varname] = e_z


def calc_lab_ehor(sensor):
    """ Calculate the horizontal unit vectors of the lab frame.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the sensor data.
    """
    e_z = sensor.data["lab_ez"]
    g_vec = np.array([
        sensor.data["iss_gx"][0],
        sensor.data["iss_gy"][0],
        sensor.data["iss_gz"][0]
        ])

    g = np.linalg.norm(g_vec)

    e_z = g_vec / g

    e_x = np.array([1, 0, 0]) - np.dot([1, 0, 0], e_z)*e_z
    # e_x = np.cross([1, 0, 0], e_z)
    e_x /= np.linalg.norm(e_x)

    sensor.data["lab_ex"] = e_x

    e_y = np.cross(e_z, e_x)
    e_y /= np.linalg.norm(e_y)

    sensor.data["lab_ey"] = e_y


def calc_trafo_iss_to_lab(sensor):
    """ Calculate the transformation matrix from iss to lab frame.

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the sensor data.
    """
    e_x = sensor.data["lab_ex"]
    e_y = sensor.data["lab_ey"]
    e_z = sensor.data["lab_ez"]

    matrix = np.array([
        e_x, e_y, e_z
    ])

    sensor.data["trafo_iss_to_lab"] = matrix


def iss_to_lab(sensor, varpattern, unit=None):
    """ Transform the iss accelerations with gravity removed to lab frame.

    Variable names must be given with a {} to be replaced by the axis name.
    E.g. for accelerations with gravity removed:
    varpattern = "a{}_gr"

    Parameters
    ----------
    sensor: wana.sensor.Sensor
        Object holding the sensor data.
    varpattern: str
        Variable pattern to transform.
    unit: str
        Physical unit of the variable.
    """
    iss_x = sensor.data["iss_" + varpattern.format("x")]
    iss_y = sensor.data["iss_" + varpattern.format("y")]
    iss_z = sensor.data["iss_" + varpattern.format("z")]
    iss_vec = np.array([iss_x, iss_y, iss_z])

    projection_matrix = sensor.data["trafo_iss_to_lab"]

    lab_vec = np.dot(projection_matrix, iss_vec)

    sensor.data["lab_" + varpattern.format("x")] = lab_vec[0]
    sensor.data["lab_" + varpattern.format("y")] = lab_vec[1]
    sensor.data["lab_" + varpattern.format("z")] = lab_vec[2]

    if unit is not None:
        sensor.units["lab_" + varpattern.format("x")] = unit
        sensor.units["lab_" + varpattern.format("y")] = unit
        sensor.units["lab_" + varpattern.format("z")] = unit

    analysis.calculate_norm(sensor, "lab_"+varpattern, unit=unit)
