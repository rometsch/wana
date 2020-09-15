""" Transformations of coordinate systems. """
import numpy as np
from scipy.spatial.transform import Rotation as R


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
