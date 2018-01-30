#ifndef FLUID_SOLVER_3D_H_
#define FLUID_SOLVER_3D_H_

#include <string>
#include <vector>
#include <fstream>
#include <climits>
#include <glm/glm.hpp>

#include "SimUtil.h"
#include "timing.h"


class FluidSolver3D {

private:
	//----------------------------------------------------------------------
	// Grid Attributes
	//----------------------------------------------------------------------

	// nx
	int m_gridWidth;
	// ny
	int m_gridHeight;
	// nz
	int m_gridDepth;
	// distance between each grid cell
	float m_dx;
	// grid of cell labels, size (nx, ny)
	SimUtil::Mat3Di m_label;

	// pressure and velocity are held in a MAC grid so that
	// p(i, j, k) = p_i_j_k
	// u(i, j, k) = u_i-1/2_j_k
	// v(i, j, k) = v_i_j-1/2_k
	// w(i, j, k) = w_i_j_k-1/2

	// grid of pressures, size (nx, ny, nz)
	SimUtil::Mat3Df m_p;
	// grid of vel x component, size (nx+1, ny, nz)
	SimUtil::Mat3Df m_u;
	// grid of vel y component, size (nx, ny+1, nz)
	SimUtil::Mat3Df m_v;
	// grid of vel z component, size (nx, ny, nz+1)
	SimUtil::Mat3Df m_w;
	// saved grid of vel x component for FLIP update, size (nx+1, ny, nz)
	SimUtil::Mat3Df m_uSaved;
	// saved grid of vel y component for FLIP update, size (nx, ny+1, nz)
	SimUtil::Mat3Df m_vSaved;
	// saved grid of vel z component for FLIP update, size (nx, ny, nz + 1)
	SimUtil::Mat3Df m_wSaved;

	//----------------------------------------------------------------------
	// Simulation Attributes
	//----------------------------------------------------------------------

	const int VEL_UNKNOWN = INT_MIN;
	// number of particles to seed in each cell at start of sim
	const int PARTICLES_PER_CELL = 8;
	// the amount of weight to give to PIC in PIC/FLIP update
	const float PIC_WEIGHT = 0.02f;
	// the maximum number of grid cells a particle should move when advected
	const int ADVECT_MAX = 1;
	// acceleration due to gravity
	const SimUtil::Vec3 GRAVITY = { 0.0f, -9.81f, 0.0f };
	// density of the fluid (kg/m^3)
	const float FLUID_DENSITY = 1000.0f;
	// surface threshold for marching squares
	const float SURFACE_THRESHOLD = 5.0f;

	// simulation time step
	float m_dt;

	//----------------------------------------------------------------------
	// Particle-related Members
	//----------------------------------------------------------------------

	// list of all particles in the simulation
	std::vector<SimUtil::Particle3D> *m_particles;

	//----------------------------------------------------------------------
	// For Timing purposes
	//----------------------------------------------------------------------
	
	//timing class object
	timing *m_timer;
	//definies if timing class is initialized
	const bool ENABLE_TIMING = true;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	// solver steps

	void seedParticles(int, std::vector<SimUtil::Particle3D>*);
	void labelGrid();
	void particlesToGrid();
	void extrapolateGridFluidData(SimUtil::Mat3Df &grid, int x, int y, int z, int extrapolationDepth);
	void saveVelocityGrids();
	void applyBodyForces();
	void applyPressure();
	void gridToParticles(float);
	void advectParticles(int);
	void cleanupParticles(float);

	// helper functions
	double trilinearHatKernel(SimUtil::Vec3);
	double hatFunction(double);
	std::vector<int> checkNeighbors(SimUtil::Mat3Di&, int[3], int[3], int[][3], int, int);
	bool hasNeighbors(SimUtil::Mat3Di&, int[3], int[3], int[][3], int, int);
	SimUtil::Vec3 interpVel(SimUtil::Mat3Df&, SimUtil::Mat3Df&, SimUtil::Mat3Df&, SimUtil::Vec3);
	void RK3(SimUtil::Particle3D*, SimUtil::Vec3, float, SimUtil::Mat3Df&, SimUtil::Mat3Df&, SimUtil::Mat3Df&);
	bool projectParticle(SimUtil::Particle3D *, float);
	std::vector<std::string> split(std::string str, std::string token);

public:
	/*
	Creates a new 2D fluid solver.
	Args:
	width - width of the grid to use
	height - height of the grid to use
	depth - depth of the grid to use
	dx - the grid cell width
	dt - the timestep to use
	*/
	FluidSolver3D(int, int, int, float, float);
	~FluidSolver3D();

	/*
	Initializes the solver by reading in and constructing initial
	grid based on the given initial geometry file. The solver will save particle
	data at each time step to the given output file.
	Args:
	initialGemoetryFile - name of the .txt file containing initial geometry
	*/
	void init(std::string);

	/*
	Steps the simulation forward dt.
	*/
	void step();
	/*
	Times the different algorithms in step()
	*/
	void stepTiming();
	/*
	Saves the lines of the isocontur representing the surface between water and
	air. Outputs in csv format where each line is represented by two points.
	each coordinate and each point is seperated by a space. Each line is one
	timestep.
	Args:
	LinesOut - pointer to file stream to use for output
	*/
	void saveParticleData(std::ofstream*);

	/*
	Returns the particles location as vectors
	*/
	std::vector<glm::vec2> particleData();
	/*
	Saves the average timing data for each sub-algorithm of the step algorithm
	*/
	void saveTimingData(std::ofstream*);
};

#endif //FLUID_SOLVER_3D_H_