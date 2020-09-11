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