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
	// distance between each grid cell
	float m_dx;
	// grid of cell labels, size (nx, ny)
	SimUtil::Mat2Di m_label;

	// pressure and velocity are held in a MAC grid so that
	// p(i, j, k) = p_i_j_k
	// u(i, j, k) = u_i-1/2_j_k
	// v(i, j, k) = v_i_j-1/2_k

	// grid of pressures, size (nx, ny)
	SimUtil::Mat2Df m_p;
	// grid of vel x component, size (nx+1, ny)
	SimUtil::Mat2Df m_u;
	// grid of vel y component, size (nx, ny+1)
	SimUtil::Mat2Df m_v;
	// saved grid of vel x component for FLIP update, size (nx+1, ny)
	SimUtil::Mat2Df m_uSaved;
	// saved grid of vel y component for FLIP update, size (nx, ny+1)
	SimUtil::Mat2Df m_vSaved;

	//----------------------------------------------------------------------
	// Simulation Attributes
	//----------------------------------------------------------------------

	const int VEL_UNKNOWN = INT_MIN;
	// number of particles to seed in each cell at start of sim
	const int PARTICLES_PER_CELL = 4;
	// the amount of weight to give to PIC in PIC/FLIP update
	const float PIC_WEIGHT = 0.02f;
	// the maximum number of grid cells a particle should move when advected
	const int ADVECT_MAX = 1;
	// acceleration due to gravity
	const SimUtil::Vec2 GRAVITY = { 0.0f, -9.81f };
	// density of the fluid (kg/m^3)
	const float FLUID_DENSITY = 1000.0f;
	// surface threshold for marching squares
	const float SURFACE_THRESHOLD = 5.0f;
	//defines if bSplines or hat function should be used for ParticleToGrid
	const bool USE_BSPLINE = false;

	// simulation time step
	float m_dt;

	//----------------------------------------------------------------------
	// Particle-related Members
	//----------------------------------------------------------------------

	// list of all particles in the simulation
	std::vector<SimUtil::Particle2D> *m_particles;

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

	void seedParticles(int, std::vector<SimUtil::Particle2D>*);
	void labelGrid();
	void particlesToGrid();
	void extrapolateGridFluidData(SimUtil::Mat2Df&, int, int, int);
	void saveVelocityGrids();
	void applyBodyForces();
	void applyPressure();
	void gridToParticles(float);
	void advectParticles(int);
	void cleanupParticles(float);

	// helper functions
	double trilinearHatKernel(SimUtil::Vec2);
	double hatFunction(double);
	double quadBSplineKernel(SimUtil::Vec2);
	double bSplineFunction(double);
	std::vector<int> checkNeighbors(SimUtil::Mat2Di&, int[2], int[2], int[][2], int, int);
	bool hasNeighbors(SimUtil::Mat2Di&, int[2], int[2], int[][2], int, int);
	SimUtil::Vec2 interpVel(SimUtil::Mat2Df&, SimUtil::Mat2Df&, SimUtil::Vec2);
	void RK3(SimUtil::Particle2D*, SimUtil::Vec2, float, SimUtil::Mat2Df&, SimUtil::Mat2Df&);
	bool projectParticle(SimUtil::Particle2D *, float);
	std::vector<glm::vec2> marchingSquares(SimUtil::Mat2Df& grid, int width, int height, float tol);
	SimUtil::MarchingTrianglesData marchingSquaresTriangles(SimUtil::Mat2Df& grid, int width, int height, float tol);
	std::vector<std::string> split(std::string str, std::string token);

public:
	/*
	Creates a new 2D fluid solver.
	Args:
	width - width of the grid to use
	height - height of the grid to use
	dx - the grid cell width
	dt - the timestep to use
	*/
	FluidSolver3D(int, int, float, float);
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
	void saveLineData(std::ofstream*);
	void saveLineData(std::ofstream*, float);
	void saveTriangleData(std::ofstream* vert, std::ofstream* ind, std::ofstream* opa);
	void saveTriangleData(std::ofstream* vert, std::ofstream* ind, std::ofstream* opa, float);

	/*
	Returns the particles location as vectors
	*/
	std::vector<glm::vec2> particleData();
	/*
	Returns the line data of the isocontour in a given grid as a vertices.
	each line has a startpoint and an end point stored as a vec2 in the vector
	*/
	std::vector<glm::vec2> marchingSquares();
	/*
	Returns the triangle data of the isocontour in the current pressure grid as a struct of 3 vectors.
	first containing vertices as vec2
	second containing vertex indices for index buffering
	third containing opacity values per vertex for pressure dependend opacity
	*/
	SimUtil::MarchingTrianglesData marchingSquaresTriangles();
	/*
	Saves the average timing data for each sub-algorithm of the step algorithm
	*/
	void saveTimingData(std::ofstream*);
};

#endif //FLUID_SOLVER_3D_H_