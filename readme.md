# wana - Walk ANAlysis

A framework for walking analysis using acceleration and gyroscope data from sensors strapped to the feet.

This framework is currently using raw data from two commercially available sensors.

The framework includes a python library which implements classes and functions to handle the following aspect:

+ `session.py`: The `Sessions` class supports loading and parsing of `session.xml` files which are writen by the proprietary sensors. The files also includes a command line part to print out a session file.
+ `task.py`: The proprietary software supports running multiple tasks/assignments in a row. The data for a single assignment is contained within the `Task` class. It holds references to the sensor objects.
+ `sensor.py`: The `Sensor` class holds the actual acceleration and angular velocity data from the sensor. It support applying a calibration file and postprocessing of the raw sensor values.
+ `calibration.py`: Interpretation of the proprietary calibration files is done in the `Calibration` class.
+ `plots.py`: Functions to plot the different values. These are used for the command line part.

## How to install

Run `./scripts/deploy.sh` on a linux maschine to install it to your user python installation or run

``` bash
python3 setup.py install --user
```

## How to use

First, make sure to copy the appropriate calibration files to the same directory as the `session.xml` file. Running `wana plot session.xml 0` should give a warning which includes the name of the necessary file.

Either import the appropriate parts from the `wana` package, e.g.,

``` python
from wana.task import Task
```

or use the command line for plotting e.g.

``` bash
wana plot session.xml 0 -p acc -o acceleration.png
```
To plot the acceleration data of the first task from the specified `session.xml` file and save it to a png file.

Run
``` bash
wana --help
wana plot --help
wana session --help
```
for further information on the usage of the command line tools.

### Plotting a specific valiable

The plot function has a manual mode which allows to visualize an arbitrary variable.
To plot all three velocities in the lab frame run

``` python
wana plot session.xml 0 -p manual -n lab_vz lab_vy lab_vz
```

## Available postprocessing data

| variable            | meaning                                                       | shape     |
| ------------------- | ------------------------------------------------------------- | --------- |
| counter             | counts the measurements done by the sensor                    | N         |
| a{x,y,z}            | accelerations in the sensor frame in m/s                      | N         |
| a                   | length of the acceleration vector                             | N         |
| r{x,y,z}            | angular velocities in deg/sec around the sensor's axis        | N         |
| time                | time in seconds                                               | N         |
| dt                  | length of timestep in seconds                                 | N         |
| delta_angle_{x,y,z} | change of angle in the interval = r{x,y,z}*dt                 | N         |
| rotation_to_iss     | rotation matrix from sensor system to iss                     | N         |
| iss_a{x,y,z}        | acceleration in the iss frame                                 | N         |
| iss_a{x,y,z}_rest   | acceleration in the iss frame but only values flagged at rest | N         |
| mask_resting        | Indicator whether sensor is resting (resting=1, moving=0)     | N         |
| mask_moving         | Indicator whether sensor is moving (resting=0, moving=1)      | N         |
| iss_g               | g vector in the iss system = average over values at rest      | 3         |
| iss_a{x,y,z}_gr     | acceleration in the iss frame with g removed                  | N         |
| iss_v{x,y,z}        | velocity in the iss frame = integral iss_a{x,y,z}_gr  dt      | N         |
| iss_{x,y,z}         | positions in the iss frame = intergral iss_v{x,y,z}_gr  dt    | N         |
| lab_e{x,y,z}        | unit vectors of the lab system in the iss frame               | 3         |
| rotation_iss_2_lab  | rotation from iss to lab frame                                | 1         |
| lab_a{x,y,z}_gr     | acceleration in lab frame with g removed                      | N         |
| lab_v{x,y,z}        | velocity in lab frame = integral lab_a{x,y,z}_gr  dt          | N         |
| lab_{x,y,z}         | position in lab frame = integral lab_v{x,y,z} dt              | N         |
| interval_steps      | a list of indices indicating start and stop of a step         | Nstep x 2 |
| iss_v{x,y,z}_step   | velocity in the iss frame, zeroed before each step            | N         |
| iss_{x,y,z}_step    | positions in the iss frame, zeroed before each step           | N         |
| lab_v{x,y,z}_step   | velocity in lab frame, zeroed before each step                | N         |
| lab_{x,y,z}_step    | position in lab frame, zeroed before each step                | N         |


## Coordinate Systems

The sensors are attached to the top of the shoes.
As a result, the coordinate system of the sensor is not aligned with the sole of the shoes or the ground.

To calculate the distance walked and equivalently the height above the ground we need to have the data in a fixed reference system. Let's call this system the lab system which is denoted by an $L$ in superscript (variables in this frame in the code have the prefix `lab`).
The original sensor system is denoted by a $S$ in superscript (code variables have no prefix).

To do this, we first need to know the initial orientation in space.
We can use the earth's acceleration as a reference when the feet are at rest.

Let us call the acceleration vector due to earth's acceleration $\vec{g}$.
Furthermore, let the unit vector in this direction be $\hat{e}_g = \frac{\vec{g}}{g}$ with $g = |\vec{g}|$.

To get the z-value in the inertial system, we simply use a projection onto the unit vector;

$$ a_z^L = \hat{e}_g^S \vec{a}^S $$

This should obviously result in the value of $g = 9.81 \frac{\mathrm{m}}{\mathrm{s}^2}$.

Since the feet are free to rotate freely, the direction of $\vec{g}$ in the sensor system is constantly changing.
Using the gyroscope data which provides angular velocities, we can compute the instantaneous $\vec{g}_n$ for any timestep $n$ by successively correcting for the rotation.

To track the orientation of the gyroscope sensor system we can use the simultaneous angular rotation angle (SORA) method as described in [Stancin & Tomazic 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3231471/).

The appropriate rotation matrix for rotating the coordinate systems is a rotation by angle $\Phi$ around a rotation vector $\vec{v}$.
Having the angular velocities $\omega^n_i$, $i \in {1,2,3}$ at a given measurement interval $n$ with duration $\Delta t^n$ the angles by which the axis are rotated are $\phi^n_i = \Delta t^n \omega^n_i$.

Assuming that the axis of rotation is constant over the measurement interval we can use the SORA method to obtain

$$ \Phi^n = \sqrt{\sum_{k=1}^3 {\phi^n_k}^2} = \Delta t^n \sqrt{\sum_{k=1}^3 {\omega^n_k}^2} $$

and 

$$ \vec{v} = \frac{1}{\Phi^n} \begin{pmatrix} \phi^n_1 \\ \phi^n_2 \\ \phi^n_3 \end{pmatrix} $$

or to express everything in a single vector

$$ \vec{\Phi^n} = \begin{pmatrix} \phi^n_1 \\ \phi^n_2 \\ \phi^n_3 \end{pmatrix}\,. $$

The length of which gives the rotation angle and its direction specifies the rotation axis.

The angles or angular velocities need to be specified in radians and radians per seconds.

To perform the 3D rotations the `scipy.transform` module ([documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html)) can be used.

``` python
from scipy.spatial.transform import Rotation as R
r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))
r.as_matrix()
```

This rotation can be computed for every timestep.
By successively applying all $n-1$ prior rotations in reverse for time-step $n$, we can calculate the rotation which transforms values from time-step $n$ to the initial sensor system at rest for timestep $n$.
We call this frame of reference the initial sensor system (ISS). Variables in the code have the `iss` prefix.

### Coordinate system conventions

| frame                     | description                                          | prefix |
| ------------------------- | ---------------------------------------------------- | ------ |
| sensor                    | values in the intrinsic sensor system                |        |
| initial sensor system     | inertial system of the sensor being at rest          | iss    |
| laboratory inertial frame | the reference system with g orthogonal to the ground | lab    |
