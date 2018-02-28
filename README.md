# 3D Real-Time Fluid Simulation

![title picture](https://raw.githubusercontent.com/a1ex90/Fluids3D/master/pictures/title.png)

This project contains a 3D grid-based real-time fluid simulation with marker particles in a hybrid pic-flip blend with a particle and a marching cubes fluid visualization.

For more detail take a look at my [paper](https://github.com/a1ex90/Fluids3D/blob/master/documentation/Sommer%2C3D-RT-Fluid-Sim.pdf).

## Install

The following include and library dependencies need to be included in the project:

* [glew](http://glew.sourceforge.net/) - OpenGL Extension Wrangler Library
* [glm](https://glm.g-truc.net/) - OpenGL Mathematics
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
| m      | De-/Active Gravity Orientation Update |
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

A surface-threshold constant in Pa (N*m) used for interpolation in marching cubes

```
const float SURFACE_THRESHOLD = ...
```

## Initial Geometry File Format

As mentioned above you're able to use own geometries containing solid-, fluid- or empty grid points. The file has to formated in a special way explained in the following example:

```
3  //number of solid grid points around the scene as a border
3  //number of empty z-layers (just containing air) in the front of the scene
3  //number of empty z-layers in the back of the scene
16 //number of scene repetitions in z-layers
//next line two dimensional scene that gets repeated f=fluid, a=empty(air), s=solid
aaaaaaaaaaaaaaaaaaa
//... as many lines as you want your grid to be high
fffffffffaaaaaaaaaa 
```

This example generates a 25x25x25 grid. Please leave out all the comments in your file.

## Authors

* **Alex Sommer** - *Initial work*

## License

This project is licensed under the MIT License.

## Acknowledgments

* Using methods found in Robert Bridson's "Fluid Simulation for Computer Graphics"
