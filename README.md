# Real-Time 3D Fluid Simulation

![title picture](https://raw.githubusercontent.com/a1ex90/Fluids3D/master/pictures/title.png)

This project contains a 3D grid-based real-time fluid simulation with marker particles in a hybrid pic-flip blend with a particle and a marching cubes fluid visualization.

## Install

The following include and library dependencies need to be included in the project:

* [glew](http://glew.sourceforge.net/) - OpenGL Extension Wrangler Library
* [glm](https://maven.apache.org/) - OpenGL Mathematics
* [SDL](https://www.libsdl.org/) - Simple DirectMedia Layer

## Usage

In the live-demo you're able to use the following input commands:

| Hotkey | Action                    |
| ------ | ------------------------- |
| e      | Zoom In                   |
| q      | Zoom Out                  |
| w      | Move Screen Up            |
| s      | Move Screen Down          |
| a      | Move Screen Left          |
| d      | Move Screen Right         |
| space  | Pause/Resume Simulation   |
| n      | Jump to Next Sim Step     |
| m      | De-/Active Gravity Rotate |
| o      | + Rotation Z-Axis         |
| l      | - Rotation Z-Axis         |
| right  | + Rotation Y-Axis         |
| left   | - Rotation Y-Axis         |
| up     | + Rotation X-Axis         |
| down   | - Rotation X-Axis         |

Furthermore you're able to do an arcball rotation with mouse dragging.

## Configurations

You're able to change the following parameters to tweak the simulation behavior in code:

### Inside FluidSim3dMain.cpp

Size of the grid-cells in meters

```
const float GRID_CELL_WIDTH = ...
```

Simulation time step, which is also the duration each frame is displayed

```
const float TIME_STEP = ...
```

Initial geometry file, explaination about the file-format later

```
const std::string INITIAL_GEOMETRY_FILE_IN = ...
```

### Inside FluidSolver3d.h

Number of particles that will initially spawn in each fluid cell 

```
const int PARTICLES_PER_CELL = ...
```

Acceleration due to gravitation in m/s^2

```
const float GRAVITY = ...
```

Density of the fluid in kg/m^3

```
const float FLUID_DENSITY = ...
```

A surface-threshold constant in N*m used for interpolation in marching cubes

```
const float SURFACE_THRESHOLD = ...
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
