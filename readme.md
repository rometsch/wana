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

## Coordinate Systems

The sensors are attached to the top of the shoes.
As a result, the coordinate system of the sensor is not aligned with the sole of the shoes or the ground.

To calculate the distance walked and equivalently the height above the ground we need to have the data in a fixed reference system. Let's call this system the inertial system which is denoted by an $I$ in superscript.
The original sensor system is denoted by a $S$ in superscript.

To do this, we first need to know the initial orientation in space.
We can use the earth's acceleration as a reference when the feet are at rest.

Let us call the acceleration vector due to earth's acceleration $\vec{g}$.
Furthermore, let the unit vector in this direction be $\hat{e}_g = \frac{\vec{g}}{g}$ with $g = |\vec{g}|$.

To get the z-value in the inertial system, we simply use a projection onto the unit vector;

$$ a_z^I = \hat{e}_g^S \vec{a}^S $$

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
