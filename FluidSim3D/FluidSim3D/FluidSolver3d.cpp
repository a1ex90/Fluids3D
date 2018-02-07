#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
//Timing
#include <chrono>
#include <ratio>

#include "FluidSolver3d.h"
#include "SimUtil.h"
#include "classicSolver.h"

using namespace SimUtil;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

FluidSolver3D::FluidSolver3D(int width, int height, int depth, float dx, float dt){
	m_gridWidth = width;
	m_gridHeight = height;
	m_gridDepth = depth;
	if (width > height) {
		if (width > depth)
			m_maxGridSize = width;
		else
			m_maxGridSize = depth;
	}
	else {
		if (height > depth)
			m_maxGridSize = height;
		else
			m_maxGridSize = depth;
	}
	m_dx = dx;
	m_dt = dt;

	m_particles = new std::vector<Particle3D>();
	if (ENABLE_TIMING) {
		std::vector<std::string> algorithms{
			"labelGrid",
			"particlesToGrid",
			"extrapolateDepth2",
			"saveVelocityGrids",
			"applyBodyForces",
			"pressureSolve",
			"applyPressure",
			"extrapolateDepthN",
			"gridToParticles",
			"advectParticles",
			"cleanupParticles"
		};
		m_timer = new timing{ algorithms };
	}
}

//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------

FluidSolver3D::~FluidSolver3D(){
	m_label.deleteGrid();
	m_p.deleteGrid();
	m_u.deleteGrid();
	m_uSaved.deleteGrid();
	m_v.deleteGrid();
	m_vSaved.deleteGrid();
	m_w.deleteGrid();
	m_wSaved.deleteGrid();

	delete m_particles;
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidSolver3D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	m_label = Mat3Di(m_gridWidth, m_gridHeight, m_gridDepth);
	m_p = Mat3Df(m_gridWidth, m_gridHeight, m_gridDepth);
	m_u = Mat3Df(m_gridWidth + 1, m_gridHeight, m_gridDepth);
	m_uSaved = Mat3Df(m_gridWidth + 1, m_gridHeight, m_gridDepth);
	m_v = Mat3Df(m_gridWidth, m_gridHeight + 1, m_gridDepth);
	m_vSaved = Mat3Df(m_gridWidth, m_gridHeight + 1, m_gridDepth);
	m_w = Mat3Df(m_gridWidth, m_gridHeight, m_gridDepth + 1);
	m_wSaved = Mat3Df(m_gridWidth, m_gridHeight, m_gridDepth + 1);

	// init vel grids with unknown label value
	m_u.initValues(VEL_UNKNOWN);
	m_v.initValues(VEL_UNKNOWN);
	m_w.initValues(VEL_UNKNOWN);

	// read in initial geometry to populate label grid
	readInGeom3D(m_gridWidth, m_gridHeight, m_gridDepth, initialGeometryFile, m_label);

	// seed particles using label grid
	seedParticles(PARTICLES_PER_CELL, m_particles);
}

void FluidSolver3D::step() {
	// update the grid labels
	labelGrid();
	// transfer particle vel to grid
	particlesToGrid();
	// extrapolate fluid data out one cell for accurate divergence calculations
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridDepth, 2);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridDepth, 2);
	extrapolateGridFluidData(m_w, m_gridWidth, m_gridHeight, m_gridDepth + 1, 2);
	// save copy of current grid velocities for FLIP update
	saveVelocityGrids();
	// apply body forces on grid (gravity)
	applyBodyForces();
	// solve for pressure
	classicSolver solver(m_gridWidth, m_gridHeight, m_gridDepth, m_dx, m_dt, &m_label, &m_p, &m_u, &m_v, &m_w);
	solver.pressureSolve();
	// apply pressure force
	applyPressure();	
	// extrapolate fluid data for the whole grid
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridDepth, m_maxGridSize);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridDepth, m_maxGridSize);
	extrapolateGridFluidData(m_w, m_gridWidth, m_gridHeight, m_gridDepth + 1, m_maxGridSize);
	// transfer grid velocities back to particles
	gridToParticles(PIC_WEIGHT);
	// advect particles
	advectParticles(ADVECT_MAX);
	// detect particles that have penetrated solid boundary and move back inside fluid
	cleanupParticles(m_dx / 4.0f);

}

void FluidSolver3D::stepTiming() {
	// update the grid labels
	m_timer->start();
	labelGrid();
	m_timer->stop();
	// transfer particle vel to grid
	m_timer->start();
	particlesToGrid();
	m_timer->stop();
	// extrapolate fluid data out one cell for accurate divergence calculations
	m_timer->start();
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridDepth, 2);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridDepth, 2);
	extrapolateGridFluidData(m_w, m_gridWidth, m_gridHeight, m_gridDepth + 1, 2);
	m_timer->stop();
	// save copy of current grid velocities for FLIP update
	m_timer->start();
	saveVelocityGrids();
	m_timer->stop();
	// apply body forces on grid (gravity)
	m_timer->start();
	applyBodyForces();
	m_timer->stop();
	// solve for pressure
	m_timer->start();
	classicSolver solver(m_gridWidth, m_gridHeight, m_gridDepth, m_dx, m_dt, &m_label, &m_p, &m_u, &m_v, &m_w);
	solver.pressureSolve();
	m_timer->stop();
	// apply pressure force
	m_timer->start();
	applyPressure();
	m_timer->stop();
	// advect particles
	m_timer->start();
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridDepth, m_maxGridSize);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridDepth, m_maxGridSize);
	extrapolateGridFluidData(m_w, m_gridWidth, m_gridHeight, m_gridDepth + 1, m_maxGridSize);
	m_timer->stop();
	// transfer grid velocities back to particles
	m_timer->start();
	gridToParticles(PIC_WEIGHT);
	m_timer->stop();
	// advect particles
	m_timer->start();
	advectParticles(ADVECT_MAX);
	m_timer->stop();
	// detect particles that have penetrated solid boundary and move back inside fluid
	m_timer->start();
	cleanupParticles(m_dx / 4.0f);
	m_timer->stop();

	m_timer->endStep();
}

void FluidSolver3D::saveParticleData(std::ofstream *particleOut) {
	if (particleOut->is_open()) {
				// print out all particle data on same line, each pos separated by " "
			size_t numParticles = m_particles->size();
		if (numParticles > 0) {
			for (int i = 0; i < numParticles - 1; i++) {
				(*particleOut) << 2 * m_particles->at(i).pos.x / (m_gridWidth * m_dx) - 1 << " " << 2 * m_particles->at(i).pos.y / (m_gridHeight * m_dx) - 1 << " ";				
			}
			(*particleOut) << 2 * m_particles->at(numParticles - 1).pos.x / (m_gridWidth * m_dx) - 1 << " " << 2 * m_particles->at(numParticles - 1).pos.y / (m_gridHeight * m_dx) - 1 << "\n";
			
		}
		else {
			(*particleOut) << "\n";
			
		}		
	}	
}

void FluidSolver3D::saveTimingData(std::ofstream *timingOut) {
	m_timer->writeTiming(timingOut);
}

std::vector<glm::vec3> FluidSolver3D::particleData() {
	std::vector<glm::vec3> particles;
	particles.reserve(m_particles->size());
	for (int i = 0; i < m_particles->size(); i++) {
		float x = 2 * m_particles->at(i).pos.x / (m_maxGridSize * m_dx) - 1;
		float y = 2 * m_particles->at(i).pos.y / (m_maxGridSize * m_dx) - 1;
		float z = 2 * m_particles->at(i).pos.z / (m_maxGridSize * m_dx) - 1;
		glm::vec3 p_i{ x,y,z };
		particles.push_back(p_i);
	}
	return particles;
}

Mesh3D FluidSolver3D::meshData() {
	return meshData(m_p, m_cubeCases, m_cubeIndices, m_gridWidth, m_gridHeight, m_gridDepth, SURFACE_THRESHOLD);
}


//----------------------------------------------------------------------
// Private Main Solver Step Functions
//----------------------------------------------------------------------

/*
Seeds the initial simulation particles. Particles are created for each fluid-labeled
cell in a random-jittered subgrid pattern.
Args:
particlesPerCell - number of particles to seed in each fluid cell.
particleList - list to place the new particles in
*/
void FluidSolver3D::seedParticles(int particlesPerCell, std::vector<Particle3D> *particleList) {
	// set random seed
	srand(time(NULL));
	// go through all cells marked fluid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label.get(i, j, k) == SimUtil::FLUID) {
					// seed randomly in 2x2 subgrid of the cell
					Vec3 cellCenter = getGridCellPosition(i, j, k, m_dx);
					Vec3 subCenters[] = {
						Vec3(cellCenter.x - 0.25f*m_dx, cellCenter.y + 0.25f*m_dx, cellCenter.z -0.25f*m_dx), // top left
						Vec3(cellCenter.x + 0.25f*m_dx, cellCenter.y + 0.25f*m_dx, cellCenter.z - 0.25f*m_dx), // top right
						Vec3(cellCenter.x + 0.25f*m_dx, cellCenter.y - 0.25f*m_dx, cellCenter.z - 0.25f*m_dx), // bottom right
						Vec3(cellCenter.x - 0.25f*m_dx, cellCenter.y - 0.25f*m_dx, cellCenter.z - 0.25f*m_dx), // bottom left
						Vec3(cellCenter.x - 0.25f*m_dx, cellCenter.y + 0.25f*m_dx, cellCenter.z + 0.25f*m_dx), // top left
						Vec3(cellCenter.x + 0.25f*m_dx, cellCenter.y + 0.25f*m_dx, cellCenter.z + 0.25f*m_dx), // top right
						Vec3(cellCenter.x + 0.25f*m_dx, cellCenter.y - 0.25f*m_dx, cellCenter.z + 0.25f*m_dx), // bottom right
						Vec3(cellCenter.x - 0.25f*m_dx, cellCenter.y - 0.25f*m_dx, cellCenter.z + 0.25f*m_dx) // bottom left
					};
					// cycle through subgrid to place all particles
					for (int k = 0; k < particlesPerCell; k++) {
						// randomly jitter from subgrid center
						// give a random factor from [-0.24, 0.24] multiplied by dx
						float jitterX = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
						//jitterX = 0.0f; //DEBUG ALEX
						float jitterY = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
						//jitterY = 0.0f; //DEBUG ALEX
						float jitterZ = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
						//jitterZ = 0.0f; //DEBUG ALEX
						Vec3 pos(subCenters[k % 8].x + jitterX, subCenters[k % 8].y + jitterY, subCenters[k % 8].z + jitterZ);
						Vec3 vel(0.0f, 0.0f, 0.0f);
						particleList->push_back(Particle3D(pos, vel));
					}
				}
			}		
		}
	}
}

/*
Updates the label grid based on particle positions. Grid cells containing particles
are fluid, solids remain the same, and all other are air. 
*/
void FluidSolver3D::labelGrid() {
	// first clear grid labels (mark everything air, but leave solids)
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label.get(i, j, k) != SOLID) {
					m_label.set(i, j, k, AIR);
				}
			}
		}
	}

	// mark any cell containing a particle FLUID
	for (int i = 0; i < m_particles->size(); i++) {
		// get cell containing the particle
		int *cell = getGridCellIndex(m_particles->at(i).pos, m_dx);
		m_label.set(cell[0], cell[1], cell[2], FLUID);
	}
}

/*
Updates the velocity field from the velocity values of the surrounding particles (PIC)
by summing up weighted partciles values at each grid point
*/
void FluidSolver3D::particlesToGrid() {
	// init vel grids with unknown label value
	m_u.initValues(VEL_UNKNOWN);
	m_v.initValues(VEL_UNKNOWN);
	m_w.initValues(VEL_UNKNOWN);

	// For each component of velocity in each fluid grid cell
	// we calculate weighted average of particles around it defined
	// by a kernel function and set this as the vel in the cell. 

	// structures to accumulate grid numerator and denominator of weighted average before divide
	Mat3Dd uNum{ m_gridWidth + 1, m_gridHeight, m_gridDepth };
	Mat3Dd uDen{ m_gridWidth + 1, m_gridHeight, m_gridDepth };
	Mat3Dd vNum{ m_gridWidth, m_gridHeight + 1, m_gridDepth };
	Mat3Dd vDen{ m_gridWidth, m_gridHeight + 1, m_gridDepth };
	Mat3Dd wNum{ m_gridWidth, m_gridHeight, m_gridDepth + 1 };
	Mat3Dd wDen{ m_gridWidth, m_gridHeight, m_gridDepth + 1 };

	// clear accumulators TODO: could be faster by initialization in one loop
	uNum.initValues(0.0);
	uDen.initValues(0.0);
	vNum.initValues(0.0);
	vDen.initValues(0.0);
	wNum.initValues(0.0);
	wDen.initValues(0.0);

	// loop over particles and accumulate num and den at each grid point
	for (int p = 0; p < m_particles->size(); p++) {
		Particle3D curParticle = m_particles->at(p);
		int* ind = getGridCellIndex(curParticle.pos, m_dx);
		int i0 = ind[0];
		int j0 = ind[1];
		int k0 = ind[2];
		for (int i = i0 - 1; i <= i0 + 1; i++) {
			for (int j = j0 - 1; j <= j0 + 1; j++) {
				for (int k = k0 - 1; k <= k0 + 1; k++) {
					if (j < m_gridHeight && k < m_gridDepth) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i - 0.5f, j, k, m_dx)));
						uNum.set(i, j, k, (uNum.get(i, j, k) + curParticle.vel.x * kernel));
						uDen.set(i, j, k, uDen.get(i, j, k) + kernel);
					}
					if (i < m_gridWidth && k < m_gridDepth) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i, j - 0.5f, k, m_dx)));
						vNum.set(i, j, k, (vNum.get(i, j, k) + curParticle.vel.y * kernel));
						vDen.set(i, j, k, vDen.get(i, j, k) + kernel);
					}
					if (i < m_gridWidth && j < m_gridHeight) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i, j, k - 0.5f, m_dx)));
						wNum.set(i, j, k, (wNum.get(i, j, k) + curParticle.vel.z * kernel));
						wDen.set(i, j, k, wDen.get(i, j, k) + kernel);
					}
				}
			}
		}
	}

	// additional pass over grid to divide and update actual velocities
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			for (int k = 0; k < m_gridDepth + 1; k++) {
				if (j < m_gridHeight && k < m_gridDepth) {
					if (uDen.get(i, j, k) != 0.0) {
						double val = uNum.get(i, j, k) / uDen.get(i, j, k);
						m_u.set(i, j, k, val);
					}
				}
				if (i < m_gridWidth && k < m_gridDepth) {
					if (vDen.get(i, j, k) != 0.0) {
						double val = vNum.get(i, j, k) / vDen.get(i, j, k);
						m_v.set(i, j, k, val);
					}
				}
				if (i < m_gridWidth && j < m_gridHeight) {
					if (wDen.get(i, j, k) != 0.0) {
						double val = wNum.get(i, j, k) / wDen.get(i, j, k);
						m_w.set(i, j, k, val);
					}
				}
			}
		}
	}

	uNum.deleteGrid();
	uDen.deleteGrid();
	vNum.deleteGrid();
	vDen.deleteGrid();
	wNum.deleteGrid();
	wDen.deleteGrid();
}

/*
Extrapolates the data in fluid cells of the given grid out using a breadth-first
search technique.
Args:
grid - the grid with data to extrapolate
x, y, z - the grid dimensions
depth - the number of cells away from fluid cells to extrapolate to.
*/
void FluidSolver3D::extrapolateGridFluidData(Mat3Df &grid, int x, int y, int z, int depth) {
	// initialize marker array 
	Mat3Di d{ x, y, z };
	// set d to 0 for known values, max int for unknown
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				if (grid.get(i, j, k) != VEL_UNKNOWN) {
					d.set(i, j, k, 0);
				}
				else {
					d.set(i, j, k, INT_MAX);
				}
			}
		}
	}
	// define neighbors
	int numNeighbors = 26;
	int neighbors[26][3] = {
		{ -1, 1, -1 }, // top left
		{ -1, 0, -1 }, // middle left
		{ -1, -1, -1 }, // bottom left
		{ 0, 1, -1 }, // top middle
		{ 0, 0, -1 }, // middle
		{ 0, -1, -1 }, // bottom middle
		{ 1, 1, -1 }, // top right
		{ 1, 0, -1 }, // middle right
		{ 1, -1, -1 }, // bottom right

		{ -1, 1, 0 }, // top left
		{ -1, 0, 0 }, // middle left
		{ -1, -1, 0 }, // bottom left
		{ 0, 1, 0 }, // top middle
		{ 0, -1, 0 }, // bottom middle
		{ 1, 1, 0 }, // top right
		{ 1, 0, 0 }, // middle right
		{ 1, -1, 0 }, // bottom right

		{ -1, 1, 1 }, // top left
		{ -1, 0, 1 }, // middle left
		{ -1, -1, 1 }, // bottom left
		{ 0, 1, 1 }, // top middle
		{ 0, 0, 1 }, // middle
		{ 0, -1, 1 }, // bottom middle
		{ 1, 1, 1 }, // top right
		{ 1, 0, 1 }, // middle right
		{ 1, -1, 1 } // bottom right
	};
	// initialize first wavefront
	std::vector<Vec3> W;
	int dim[3] = { x, y, z };
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				// current value is not known
				if (d.get(i, j, k) != 0) {
					int ind[3] = { i, j, k };
					if (hasNeighbors(d, dim, ind, neighbors, numNeighbors, 0)) {
						// a neighbor is known
						d.set(i, j, k, 1);
						W.push_back(Vec3(i, j, k));
					}
				}
			}
		}
	}
	// list of all wavefronts, only want to go through the given depth
	std::vector<std::vector<Vec3>> wavefronts;
	wavefronts.push_back(W);
	int curWave = 0;
	while (curWave < depth) {
		// get wavefront
		std::vector<Vec3> curW = wavefronts.at(curWave);
		// initialize next wavefront
		std::vector<Vec3> nextW;
		// go through current wave and extrapolate values
		for (int i = 0; i < curW.size(); i++) {
			Vec3 ind = curW.at(i);
			// average neighbors
			float avg = 0.0f;
			int numUsed = 0;
			for (int j = 0; j < numNeighbors; j++) {
				int offsetX = neighbors[j][0];
				int offsetY = neighbors[j][1];
				int offsetZ = neighbors[j][2];
				int neighborX = ind.x + offsetX;
				int neighborY = ind.y + offsetY;
				int neighborZ = ind.z + offsetZ;

				// make sure valid indices
				if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1]) && (neighborZ >= 0 && neighborZ < dim[2])) {
					// only want to add to average if neighbor d is less than current d
					if (d.get(neighborX, neighborY, neighborZ) < d.get((int)ind.x, (int)ind.y, (int)ind.z)) {
						avg += grid.get(neighborX, neighborY, neighborZ);
						numUsed++;
					} else if (d.get(neighborX, neighborY, neighborZ) == INT_MAX) {
						d.set(neighborX, neighborY, neighborZ, d.get((int)ind.x, (int)ind.y, (int)ind.z) + 1);
						nextW.push_back(Vec3(neighborX, neighborY, neighborZ));
					}
				}
			}

			avg /= numUsed;
			// set current value to average of neighbors
			grid.set((int)ind.x, (int)ind.y, (int)ind.z, avg);
		}

		// push next wave to list
		wavefronts.push_back(nextW);
		curWave++;
	}
	
	d.deleteGrid();
}

/*
Saves a copy of the current velocity grids to be used on the FLIP updated
*/
void FluidSolver3D::saveVelocityGrids() {
	// save u grid
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				m_uSaved.set(i, j, k, m_u.get(i, j, k));
			}
		}
	}

	// save v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				m_vSaved.set(i, j, k, m_v.get(i, j, k));
			}
		}
	}

	// save w grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth + 1; k++) {
				m_wSaved.set(i, j, k, m_w.get(i, j, k));
			}
		}
	}
}

/*
Applies the force of gravity to velocity field on the grid
*/
void FluidSolver3D::applyBodyForces() {
	// traverse all grid cells and apply force to each velocity component
	// The new velocity is calculated using forward euler
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			for (int k = 0; k < m_gridDepth + 1; k++) {
				if (j < m_gridHeight && k < m_gridDepth) {
					// make sure we know the velocity
					if (m_u.get(i, j, k) != VEL_UNKNOWN) {
						// update u component
						m_u.set(i, j, k, m_u.get(i, j, k) + m_dt*GRAVITY.x);
					}
				}
				if (i < m_gridWidth && k < m_gridDepth) {
					if (m_v.get(i, j, k) != VEL_UNKNOWN) {
						// update v component
						m_v.set(i, j, k, m_v.get(i, j, k) + m_dt*GRAVITY.y);
					}
				}
				if (i < m_gridWidth && j < m_gridHeight) {
					if (m_w.get(i, j, k) != VEL_UNKNOWN) {
						// update v component
						m_w.set(i, j, k, m_w.get(i, j, k) + m_dt*GRAVITY.z);
					}
				}
			}
		}
	}
}


/*
Applies the pressure force to the current velocity field.
*/
void FluidSolver3D::applyPressure() {
	float scale = m_dt / (FLUID_DENSITY * m_dx);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				// update u
				if (i - 1 >= 0) {
					if (m_label.get(i - 1, j, k) == FLUID || m_label.get(i, j, k) == FLUID) {
						if (m_label.get(i - 1, j, k) == SOLID || m_label.get(i, j, k) == SOLID) {
							m_u.set(i, j, k, 0.0f); // usolid[i][j]
						}
						else {
							m_u.set(i, j, k, m_u.get(i, j, k) - scale * (m_p.get(i, j, k) - m_p.get(i - 1, j, k)));
						}
					}
					else {
						m_u.set(i, j, k, VEL_UNKNOWN);
					}
				}
				else {
					// edge of grid, keep the same velocity
				}

				// update v
				if (j - 1 >= 0) {
					if (m_label.get(i, j - 1, k) == FLUID || m_label.get(i, j, k) == FLUID) {
						if (m_label.get(i, j - 1, k) == SOLID || m_label.get(i, j, k) == SOLID) {
							m_v.set(i, j, k, 0.0f); // vsolid[i][j]
						}
						else {
							m_v.set(i, j, k, m_v.get(i, j, k) - scale * (m_p.get(i, j, k) - m_p.get(i, j - 1, k)));
						}
					}
					else {
						m_v.set(i, j, k, VEL_UNKNOWN);
					}
				}
				else {
					// edge of grid, keep the same velocity
				}

				// update w
				if (k - 1 >= 0) {
					if (m_label.get(i, j, k - 1) == FLUID || m_label.get(i, j, k) == FLUID) {
						if (m_label.get(i, j, k - 1) == SOLID || m_label.get(i, j, k) == SOLID) {
							m_w.set(i, j, k, 0.0f); // wsolid[i][j]
						}
						else {
							m_w.set(i, j, k, m_v.get(i, j, k) - scale * (m_p.get(i, j, k) - m_p.get(i, j, k - 1)));
						}
					}
					else {
						m_w.set(i, j, k, VEL_UNKNOWN);
					}
				}
				else {
					// edge of grid, keep the same velocity
				}
			}
		}
	}
}

/*
Transfer the velocities from the grid back to the particles. This is done
with a PIC/FLIP mix, where the PIC update has a weight of the given alpha.
Args:
alpha - weight in the update for PIC, should be in [0, 1]. Then the FLIP update is weighted (1 - alpha).
*/
void FluidSolver3D::gridToParticles(float alpha) {
	// build grid for change in velocity to use for FLIP update
	Mat3Df duGrid{ m_gridWidth + 1, m_gridHeight, m_gridDepth };
	Mat3Df dvGrid{ m_gridWidth, m_gridHeight + 1, m_gridDepth };
	Mat3Df dwGrid{ m_gridWidth, m_gridHeight, m_gridDepth + 1 };
	// calc u grid
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				duGrid.set(i, j, k, m_u.get(i, j, k) - m_uSaved.get(i, j, k));
			}
		}
	}
	// calc v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				dvGrid.set(i, j, k, m_v.get(i, j, k) - m_vSaved.get(i, j, k));
			}
		}
	}
	// calc w grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth + 1; k++) {
				dwGrid.set(i, j, k, m_w.get(i, j, k) - m_wSaved.get(i, j, k));
			}
		}
	}

	// go through particles and interpolate each velocity component
	// the update is a PIC/FLIP mix weighted with alpha
	// alpha = 1.0 is entirely PIC, alpha = 0.0 is all FLIP
	for (int i = 0; i < m_particles->size(); i++) {
		Particle3D *curParticle = &(m_particles->at(i));
		Vec3 picInterp = interpVel(m_u, m_v, m_w, curParticle->pos);
		Vec3 flipInterp = interpVel(duGrid, dvGrid, dwGrid, curParticle->pos);
		// u_new = alpha * interp(u_gridNew, x_p) + (1 - alpha) * (u_pOld + interp(u_dGrid, x_p))
		curParticle->vel = add(scale(picInterp, alpha), scale(add(curParticle->vel, flipInterp), 1.0f - alpha));
	}

	duGrid.deleteGrid();
	dvGrid.deleteGrid();
	dwGrid.deleteGrid();
}

/*
Advects the particles through the current velocity field using a Runge-Kutta 3 method.
This uses substeps so that the particles are never advected greater than C*dx in a single
substep.
Args:
C - the maximum number of grid cells a particle should move when advected. This helps define substep sizes. 
*/
void FluidSolver3D::advectParticles(int C) {
	for (int i = 0; i < m_particles->size(); i++) {
		Particle3D *curParticle = &(m_particles->at(i));
		float subTime = 0;
		bool finished = false;
		//float dT = m_dt / 4.999f;
		while (!finished) {
			Vec3 curVel = interpVel(m_u, m_v, m_w, curParticle->pos);

			// calc max substep size
			
			float dT = (C * m_dx) / (norm(curVel) + FLT_MIN);
			// SET TEMPORARY
			//float dT = m_dt / 4.999f;
			// update substep time so we don't go past normal time step
			if (subTime + dT >= m_dt) {
				dT = m_dt - subTime;
				finished = true;
			} else if (subTime + 2 * dT >= m_dt) {
				dT = 0.5f * (m_dt - subTime);
			}

			RK3(curParticle, curVel, dT, m_u, m_v, m_w);
			subTime += dT;

			if (curParticle->pos.x < 0 || curParticle->pos.y < 0 || curParticle->pos.z < 0 ||
				 curParticle->pos.x > m_gridWidth * m_dx || curParticle->pos.y > m_gridHeight * m_dx || curParticle->pos.z > m_gridDepth * m_dx ||
				isnan(curParticle->pos.x) || isnan(curParticle->pos.y) || isnan(curParticle->pos.z)) {
				// there's been an error in RK3, just skip it
				std::cout << "RK3 error...skipping particle outside of grid" << std::endl;
				break;
			}

			int *cell = getGridCellIndex(curParticle->pos, m_dx);
			int j = cell[0];
			int k = cell[1];
			int l = cell[2];
			if (m_label.get(j, k, l) == SOLID) {
				//std::cout << "Advected into SOLID, projecting back!\n";
				if (!projectParticle(curParticle, m_dx / 4.0f)) {
					std::cout << "RK3 error...skipping particle in solid can't bring back" << std::endl;
					break;
				}

			}

		}
	}
}

/*
Finds particles that have traveled into solid boundaries and projects them back into the fluid region. 
They will be put at the closest boundary + the given dx into the fluid.
Args
dx - the amount to project stray particles away from the wall.
*/
void FluidSolver3D::cleanupParticles(float dx) {
	int i = 0;
	bool finished = false;
	int numDeleted = 0;
	while(!finished && m_particles->size() > 0) {
		int *cell = getGridCellIndex(m_particles->at(i).pos, m_dx);
		int ind[3] = { cell[0], cell[1], cell[2] };
		// if either of cells are negative or greater than sim dimensions it has left sim area
		if (ind[0] < 0 || ind[1] < 0 || ind[2] < 0 || ind[0] >= m_gridWidth || ind[1] >= m_gridHeight || ind[2] >= m_gridDepth 
			|| isnan(m_particles->at(i).pos.x) || isnan(m_particles->at(i).pos.y) || isnan(m_particles->at(i).pos.z)) {
			m_particles->erase(m_particles->begin() + i);
			numDeleted++;
			if (i >= m_particles->size()) {
				finished = true;
			}
		} else if (m_label.get(ind[0], ind[1], ind[2]) == SOLID) {
			// project back into fluid
			bool success = projectParticle(&(m_particles->at(i)), dx);
			if (!success) {
				// no near fluid, just delete
				m_particles->erase(m_particles->begin() + i);
				numDeleted++;
				if (i >= m_particles->size()) {
					finished = true;
				}
			}
			else {
				i++;
			}
		} else {
			i++;
			if (i >= m_particles->size()) {
				finished = true;
			}
		}
	}
}



//----------------------------------------------------------------------
// Private Helper Functions
//----------------------------------------------------------------------

/*
Checks neighbors of the given index in the given grid for the given value. Returns a vector
of neighbor indices (row index in the given neighbor array) that have the given value.
Args
grid - the 2D grid to look in
dim - the grid dimensions [x y]
index - the (i, j) index of the cell to look around
neighbors - the definintion of neighbors, this is a n x 2 array where each row is a pair of offsets from the index of interest
numNeighbors - the number of neighbors
value - the value to look for
*/
std::vector<int> FluidSolver3D::checkNeighbors(Mat3Di &grid, int dim[3], int index[3], int neighbors[][3], int numNeighbors, int value) {
	std::vector<int> neighborsTrue;
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int offsetZ = neighbors[i][2];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;
		int neighborZ = index[2] + offsetZ;

		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1]) && (neighborZ >= 0 && neighborZ < dim[2])) {
			if (grid.get(neighborX, neighborY, neighborZ) == value) {
				neighborsTrue.push_back(i);
			}
		}
	}

	return neighborsTrue;
}

/*
Checks neighbors of the given index in the given grid for the given value. Returns true if
neighbors with given value exists, false otherwise
Args
grid - the 2D grid to look in
dim - the grid dimensions [x y]
index - the (i, j) index of the cell to look around
neighbors - the definintion of neighbors, this is a n x 2 array where each row is a pair of offsets from the index of interest
numNeighbors - the number of neighbors
value - the value to look for
*/
bool FluidSolver3D::hasNeighbors(Mat3Di &grid, int dim[3], int index[3], int neighbors[][3], int numNeighbors, int value) {
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int offsetZ = neighbors[i][2];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;
		int neighborZ = index[2] + offsetZ;
		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1]) && (neighborZ >= 0 && neighborZ < dim[2])) {
			if (grid.get(neighborX, neighborY, neighborZ) == value) {
				return true;
			}
		}
	}
	return false;
}


/*
Returns the value of the trilinear hat function for the given
distance (x, y).
*/
double FluidSolver3D::trilinearHatKernel(SimUtil::Vec3 dist) {
	return hatFunction(dist.x / m_dx) * hatFunction(dist.y / m_dx) * hatFunction(dist.z / m_dx);
}

/*
Calculates the value of hat function for the given r. 
*/
double FluidSolver3D::hatFunction(double r) {
	double rAbs = abs(r);
	if (rAbs <= 1) {
		return 1.0 - rAbs;
	} else {
		return 0.0;
	}
}

/*
Interpolates the value in the given velocity grid at the given position using bilinear interpolation.
Returns velocity unkown if position is not on simulation grid. 
Args:
uGrid - the u component grid to interpolate from
vGrid - the v component grid to interpolate from
wGrid - the w component grid to interpolate from
pos - the position to interpolate at
*/
Vec3 FluidSolver3D::interpVel(SimUtil::Mat3Df &uGrid, SimUtil::Mat3Df &vGrid, SimUtil::Mat3Df &wGrid, Vec3 pos) {
	// get grid cell containing position
	int *cell = getGridCellIndex(pos, m_dx);
	int i = cell[0];
	int j = cell[1];
	int k = cell[2];
	// make sure this is a valid index
	if (i >= 0 && i < m_gridWidth && j >= 0 && j < m_gridHeight && k >= 0 && k < m_gridDepth) {
		// get positions of u, v, w component stored on each side of cell
		Vec3 cellLoc = getGridCellPosition(i, j, k, m_dx);
		float offset = m_dx / 2.0f;
		float x1 = cellLoc.x - offset;
		float x2 = cellLoc.x + offset;
		float y1 = cellLoc.y - offset;
		float y2 = cellLoc.y + offset;
		float z1 = cellLoc.z - offset;
		float z2 = cellLoc.z + offset;
		// get actual values at these positions
		float u1 = uGrid.get(i, j, k);
		float u2 = uGrid.get(i + 1, j, k);
		float v1 = vGrid.get(i, j, k);
		float v2 = vGrid.get(i, j + 1, k);
		float w1 = wGrid.get(i, j, k);
		float w2 = wGrid.get(i, j, k + 1);


		// the interpolated values
		float u = ((x2 - pos.x) / (x2 - x1)) * u1 + ((pos.x - x1) / (x2 - x1)) * u2;
		float v = ((y2 - pos.y) / (y2 - y1)) * v1 + ((pos.y - y1) / (y2 - y1)) * v2;
		float w = ((z2 - pos.z) / (z2 - z1)) * w1 + ((pos.z - z1) / (z2 - z1)) * w2;

		return Vec3(u, v, w);
	} else {
		return Vec3(VEL_UNKNOWN, VEL_UNKNOWN, VEL_UNKNOWN);
	}
}

/*
Advects a particle using Runge-Kutta 3 through the given velocity field.
Args:
particle - the particle to advect
initVel - the particles initial velocity in the current field, can leave UNKNOWN
dt - the time step
uGrid/vGrid/wGrid - the velocity grids to advect through
*/
void FluidSolver3D::RK3(SimUtil::Particle3D *particle, SimUtil::Vec3 initVel, float dt, SimUtil::Mat3Df &uGrid, SimUtil::Mat3Df &vGrid, SimUtil::Mat3Df &wGrid) {
	if (initVel.x == VEL_UNKNOWN && initVel.y == VEL_UNKNOWN && initVel.z == VEL_UNKNOWN) {
		initVel = interpVel(uGrid, vGrid, wGrid, particle->pos);
	}

	Vec3 k1 = initVel;
	Vec3 k2 = interpVel(uGrid, vGrid, wGrid, add(particle->pos, scale(k1, 0.5f*dt)));
	Vec3 k3 = interpVel(uGrid, vGrid, wGrid, add(particle->pos, scale(k2, 0.75f*dt)));
	k1 = scale(k1, (2.0f / 9.0f)*dt);
	k2 = scale(k2, (3.0f / 9.0f)*dt);
	k3 = scale(k3, (4.0f / 9.0f)*dt);

	particle->pos = add(particle->pos, add(k1, add(k2, k3)));
}

/*
Projects a particle within a solid back into the closest fluid or air. Returns true if
successful and false otherwise.
Args
particle - the particle to project. 
dx - the amount to project the particle from the solid wall
*/

bool FluidSolver3D::projectParticle(Particle3D *particle, float dx) {
	// project back into fluid
	// find neighbors that are fluid
	// define neighbors
	int numNeighbors = 26;
	int neighbors[26][3] = {
		{ -1, 1, -1 }, // top left
		{ -1, 0, -1 }, // middle left
		{ -1, -1, -1 }, // bottom left
		{ 0, 1, -1 }, // top middle
		{ 0, 0, -1 }, // middle
		{ 0, -1, -1 }, // bottom middle
		{ 1, 1, -1 }, // top right
		{ 1, 0, -1 }, // middle right
		{ 1, -1, -1 }, // bottom right

		{ -1, 1, 0 }, // top left
		{ -1, 0, 0 }, // middle left
		{ -1, -1, 0 }, // bottom left
		{ 0, 1, 0 }, // top middle
		{ 0, -1, 0 }, // bottom middle
		{ 1, 1, 0 }, // top right
		{ 1, 0, 0 }, // middle right
		{ 1, -1, 0 }, // bottom right

		{ -1, 1, 1 }, // top left
		{ -1, 0, 1 }, // middle left
		{ -1, -1, 1 }, // bottom left
		{ 0, 1, 1 }, // top middle
		{ 0, 0, 1 }, // middle
		{ 0, -1, 1 }, // bottom middle
		{ 1, 1, 1 }, // top right
		{ 1, 0, 1 }, // middle right
		{ 1, -1, 1 } // bottom right
	};
	int dim[3] = { m_gridWidth, m_gridHeight, m_gridDepth };
	int *cell = getGridCellIndex(particle->pos, m_dx);
	int index[3] = { cell[0], cell[1], cell[2] };
	// get neighbors that are fluid
	std::vector<int> neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, FLUID);
	if (neighborInd.size() == 0) {
		// try with air
		neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, AIR);
	}
	// find closest to particle
	int closestInd = -1;
	float closestDist = std::numeric_limits<float>::max();
	Vec3 closestVec(0.0f, 0.0f, 0.0f);
	for (int j = 0; j < neighborInd.size(); j++) {
		// get vec from particle to neighbor ind
		int ind[3] = { index[0] + neighbors[neighborInd.at(j)][0], index[1] + neighbors[neighborInd.at(j)][1], index[2] + neighbors[neighborInd.at(j)][2] };
		Vec3 cellPos = getGridCellPosition(ind[0], ind[1], ind[2], m_dx);
		Vec3 distVec = sub(cellPos, particle->pos);
		float dist = norm(distVec);
		if (dist < closestDist) {
			closestDist = dist;
			closestInd = neighborInd.at(j);
			closestVec = distVec;
		}
	}

	if (closestInd == -1) {
		return false;
	}
	else {
		// project different ways based on where closest neighbor is
		// also make sure to only project the amount given
		Vec3 projectVec(0.0f, 0.0f, 0.0f);
		switch (closestInd)
		{
		case 0:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 1:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 2:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 3:
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 4:
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 5:
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 6:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 7:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 8:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (dx - (m_dx / 2.0f));
			break;
		case 9:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			break;
		case 10:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			break;
		case 11:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			break;
		case 12:
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			break;
		case 13:
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			break;
		case 14:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			break;
		case 15:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			break;
		case 16:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			break;
		case 17:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 18:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 19:
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 20:
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 21:
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 22:
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 23:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 24:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		case 25:
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
			projectVec.z = closestVec.z + (-dx + (m_dx / 2.0f));
			break;
		default:
			break;
		}

		particle->pos = add(particle->pos, projectVec);

		return true;
	}
}

Mesh3D FluidSolver3D::meshData(SimUtil::Mat3Df &grid, std::vector<std::vector<glm::vec3>> &cubeCases, std::vector<std::vector<int>> &cubeIndices, int width, int height, int depth, float tol) {
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<int> globalIndices;
	int curInd = 0;

	for (int i = 0; i < width - 1; i++) {
		for (int j = 0; j < height - 1; j++) {
			for (int k = 0; k < depth - 1; k++) {
				//determine which of the 256 different cases for one square exists
				int selectCase = 0;
				if (grid.get(i, j, k + 1) > tol)
					selectCase += 128;
				if (grid.get(i + 1, j, k + 1) > tol)
					selectCase += 64;
				if (grid.get(i + 1, j + 1, k + 1) > tol)
					selectCase += 32;
				if (grid.get(i, j + 1, k + 1) > tol)
					selectCase += 16;
				if (grid.get(i, j, k) > tol)
					selectCase += 8;
				if (grid.get(i + 1, j, k) > tol)
					selectCase += 4;
				if (grid.get(i + 1, j + 1, k) > tol)
					selectCase += 2;
				if (grid.get(i, j + 1, k) > tol)
					selectCase += 1;
				for (int l = 0; l < cubeIndices[selectCase].size(); l++) {
					globalIndices.push_back(cubeIndices[selectCase][l] + curInd);
				}
				for (int l = 0; l < cubeCases[selectCase].size(); l++) {
					float offsetX = cubeCases[selectCase][l].x;
					float offsetY = cubeCases[selectCase][l].y;
					float offsetZ = cubeCases[selectCase][l].z;

					glm::vec3 normal;
					//the offset of 0.5f indicates that this axis needs interpolation
					if (offsetX == 0.5f) {
						int offsetJ = (int)offsetY;
						int offsetK = (int)offsetZ;
						offsetX = 1 / (grid.get(i + 1, j + offsetJ, k + offsetK) - grid.get(i, j + offsetJ, + k + offsetK)) * (tol - grid.get(i, j + offsetJ, + k + offsetK));

						normal.x = offsetX * grid.get(i, j + offsetJ, k + offsetK) - (1.0f - offsetX) * grid.get(i + 1, j + offsetJ, k + offsetK);
						normal.y = offsetX * (grid.get(i, j + offsetJ - 1, k + offsetK) - grid.get(i, j + offsetJ + 1, k + offsetK)) / 2.0 + 
							(1.0f - offsetX) * (grid.get(i + 1, j + offsetJ - 1, k + offsetK) - grid.get(i + 1, j + offsetJ + 1, k + offsetK)) / 2.0;
						normal.z = offsetX * (grid.get(i, j + offsetJ, k + offsetK - 1) - grid.get(i, j + offsetJ, k + offsetK + 1)) / 2.0 +
							(1.0f - offsetX) * (grid.get(i + 1, j + offsetJ, k + offsetK - 1) - grid.get(i + 1, j + offsetJ, k + offsetK + 1)) / 2.0;
							
					}
					else if (offsetY == 0.5f) {
						int offsetI = (int)offsetX;
						int offsetK = (int)offsetZ;
						offsetY = 1 / (grid.get(i + offsetI, j + 1, k + offsetK) - grid.get(i + offsetI, j, k + offsetK)) * (tol - grid.get(i + offsetI, j, k + offsetK));

						normal.y = offsetY * grid.get(i + offsetI, j, k + offsetK) - (1.0f - offsetY) * grid.get(i + offsetI, j + 1, k + offsetK);
						normal.x = offsetY * (grid.get(i + offsetI - 1, j, k + offsetK) - grid.get(i + offsetI + 1, j, k + offsetK)) / 2.0 +
							(1.0f - offsetY) * (grid.get(i + offsetI - 1, j + 1, k + offsetK) - grid.get(i + offsetI + 1, j + 1, k + offsetK)) / 2.0;
						normal.z = offsetY * (grid.get(i + offsetI, j, k + offsetK - 1) - grid.get(i + offsetI, j, k + offsetK + 1)) / 2.0 +
							(1.0f - offsetY) * (grid.get(i + offsetI, j + 1, k + offsetK - 1) - grid.get(i + offsetI, j + 1, k + offsetK + 1)) / 2.0;
					}
					else if (offsetZ == 0.5f) {
						int offsetI = (int)offsetX;
						int offsetJ = (int)offsetY;
						offsetY = 1 / (grid.get(i + offsetI, j + offsetJ, k + 1) - grid.get(i + offsetI, j + offsetJ, k)) * (tol - grid.get(i + offsetI, j + offsetJ, k));

						normal.z = offsetZ * grid.get(i + offsetI, j + offsetJ, k) - (1.0f - offsetZ) * grid.get(i + offsetI, j + offsetJ, k + 1);
						normal.x = offsetZ * (grid.get(i + offsetI - 1, j + offsetJ, k) - grid.get(i + offsetI + 1, j + offsetJ, k)) / 2.0 +
							(1.0f - offsetZ) * (grid.get(i + offsetI - 1, j + offsetJ, k + 1) - grid.get(i + offsetI + 1, j + offsetJ, k + 1)) / 2.0;
						normal.y = offsetZ * (grid.get(i + offsetI, j + offsetJ - 1, k) - grid.get(i + offsetI, j + offsetJ + 1, k)) / 2.0 +
							(1.0f - offsetZ) * (grid.get(i + offsetI, j + offsetJ - 1, k + 1) - grid.get(i + offsetI, j + offsetJ + 1, k + 1)) / 2.0;
					}

					float x = 2.0f * (i + offsetX) / (m_maxGridSize - 1) - 1;
					float y = 2.0f * (j + offsetY) / (m_maxGridSize - 1) - 1;
					float z = 2.0f * (j + offsetZ) / (m_maxGridSize - 1) - 1;

					//normalize normals vector
					glm::normalize(normal);

					vertices.push_back(glm::vec3{ x,y,z });
					normals.push_back(normal);

					//increment the current Index
					curInd++;
				}
			}
		}
	}
	return Mesh3D(vertices, normals, globalIndices);
}

void FluidSolver3D::initMarchingCubesCases(std::vector<std::vector<glm::vec3>> &points, std::vector<std::vector<int>> &indices) {
	points.reserve(256);
	indices.reserve(256);
	//There are 15 ambigious cases for the 256 total cases
	// ---CASE 0------------
	points[0] = std::vector<glm::vec3>{};
	indices[0] = std::vector<int>{};
	points[255] = std::vector<glm::vec3>{};
	indices[255] = std::vector<int>{};
	// ---CASE 1------------
	points[1] = std::vector<glm::vec3>{ glm::vec3 {0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
	indices[1] = std::vector<int>{ 0,1,2 };
	points[2] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[2] = std::vector<int>{ 0,1,2 };
	points[4] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[4] = std::vector<int>{ 0,1,2 };
	points[8] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[8] = std::vector<int>{ 0,1,2 };
	points[16] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[16] = std::vector<int>{ 0,1,2 };
	points[32] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[32] = std::vector<int>{ 0,1,2 };
	points[64] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[64] = std::vector<int>{ 0,1,2 };
	points[128] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[128] = std::vector<int>{ 0,1,2 };
	// ---CASE 2------------
	points[3] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[3] = std::vector<int>{ 0,1,2,1,2,3 };
	points[6] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[6] = std::vector<int>{ 0,1,2,1,2,3 };
	points[9] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[9] = std::vector<int>{ 0,1,2,1,2,3 };
	points[12] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[12] = std::vector<int>{ 0,1,2,1,2,3 };
	points[17] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[17] = std::vector<int>{ 0,1,2,1,2,3 };
	points[34] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[34] = std::vector<int>{ 0,1,2,1,2,3 };
	points[48] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[48] = std::vector<int>{ 0,1,2,1,2,3 };	
	points[68] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[68] = std::vector<int>{ 0,1,2,1,2,3 };
	points[96] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[96] = std::vector<int>{ 0,1,2,1,2,3 };
	points[136] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[136] = std::vector<int>{ 0,1,2,1,2,3 };
	points[144] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[144] = std::vector<int>{ 0,1,2,1,2,3 };
	points[192] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[192] = std::vector<int>{ 0,1,2,1,2,3 };
	// ---CASE 3------------
	points[5] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[5] = std::vector<int>{ 0,1,2,3,4,5 };
	points[10] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[10] = std::vector<int>{ 0,1,2,3,4,5 };
	points[24] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[24] = std::vector<int>{ 0,1,2,3,4,5 };
	points[38] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[38] = std::vector<int>{ 0,1,2,3,4,5 };
	points[66] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[66] = std::vector<int>{ 0,1,2,3,4,5 };
	points[80] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[80] = std::vector<int>{ 0,1,2,3,4,5 };
	points[129] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[129] = std::vector<int>{ 0,1,2,3,4,5 };
	points[160] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[160] = std::vector<int>{ 0,1,2,3,4,5 };
	// ---CASE 4------------
	points[19] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[19] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
	points[50] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[50] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

	points[140] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[140] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
	points[196] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[196] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

	points[98] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[98] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
	points[38] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[38] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

	


	// ---CASE 5------------
	points[15] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[15] = std::vector<int>{ 0,1,2,1,3,2 };
	points[51] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[51] = std::vector<int>{ 0,1,2,0,2,3 };
	points[102] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[102] = std::vector<int>{ 0,1,2,1,2,3 };
	points[153] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
	indices[153] = std::vector<int>{ 0,1,2,1,2,3 };
	points[204] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[204] = std::vector<int>{ 0,1,2,0,2,3 };
	points[240] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[240] = std::vector<int>{ 0,1,2,1,3,2 };
	
	// ---CASE 7------------
	points[90] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[90] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };
	points[165] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[165] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };


	// ---CASE 10------------
	points[20] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[20] = std::vector<int>{ 0,1,2,3,4,5 };
	points[40] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[40] = std::vector<int>{ 0,1,2,3,4,5 };
	points[65] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[65] = std::vector<int>{ 0,1,2,3,4,5 };
	points[130] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[130] = std::vector<int>{ 0,1,2,3,4,5 };
	// ---CASE 11------------
	points[21] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[21] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[22] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[22] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[28] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[28] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[41] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[41] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[42] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[42] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[44] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[44] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[52] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[52] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[56] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[56] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[67] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[67] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[69] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
	indices[69] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[73] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[73] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[81] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[81] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[84] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[84] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[97] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
	indices[97] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[104] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
	indices[104] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[131] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f },  glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[131] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[134] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[134] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[138] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[138] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[146] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[146] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[148] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
	indices[148] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[162] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[162] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[168] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[168] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[193] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
	indices[193] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	points[194] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
	indices[194] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
	// ---CASE 12------------
	points[37] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[37] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[74] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[74] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[82] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[82] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[88] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
	indices[88] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[133] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[133] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[138] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[138] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };	
	points[161] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[161] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
	points[164] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
	indices[164] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

}

/*
splits a given string at each occurring given token. returns a vector the substrings
Args:
str - input string to be splitted
token - character where the splitting should happen
*/
std::vector<std::string> FluidSolver3D::split(std::string str, std::string token) {
	std::vector<std::string>result;
	while (str.size()) {
		int index = str.find(token);
		if (index != std::string::npos) {
			result.push_back(str.substr(0, index));
			str = str.substr(index + token.size());
			if (str.size() == 0)result.push_back(str);
		}
		else {
			result.push_back(str);
			str = "";
		}
	}
	return result;
}

//----------------------------------------------------------------------
// Debugging Functions
//----------------------------------------------------------------------
void FluidSolver3D::gridValues(Mat3Df &grid, std::string name, int x, int y, int z) {
	int count = 0;
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				if (abs(grid.get(i, j, k)) > 1000) {
					//std::cout << "Strange Gridvalues at (" << i << ", " << j << ", " << k << ") of " << grid.get(i, j, k) << "\n";
					count++;
				}
			}
		}
	}
	std::cout << "Count of high velocities in "<< name << ": " << count << "\n";
}
void FluidSolver3D::strangeParticles() {
	for (int i = 0; i < m_particles->size(); i++) {
		if (abs(m_particles->at(i).pos.x) > 2 * m_gridWidth * m_dx || abs(m_particles->at(i).pos.y) > 2 * m_gridHeight * m_dx || abs(m_particles->at(i).pos.z) > 2 * m_gridDepth * m_dx) {
			std::cout << "Strange Particle NO" << i << "Pos (" << m_particles->at(i).pos.x << ", " << m_particles->at(i).pos.y << ", " << m_particles->at(i).pos.z << "\n";
		}
		if (abs(m_particles->at(i).vel.x) > 1000 || abs(m_particles->at(i).vel.y) > 1000 || abs(m_particles->at(i).vel.z) > 1000) {
			std::cout << "Strange Particle NO" << i << "Vel (" << m_particles->at(i).vel.x << ", " << m_particles->at(i).vel.y << ", " << m_particles->at(i).vel.z << "\n";
		}
	}
}

