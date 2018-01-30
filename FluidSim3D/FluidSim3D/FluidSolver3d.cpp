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

FluidSolver3D::FluidSolver3D(int width, int height, float dx, float dt){
	m_gridWidth = width;
	m_gridHeight = height;
	m_dx = dx;
	m_dt = dt;

	m_particles = new std::vector<Particle2D>();
	if (ENABLE_TIMING) {
		std::vector<std::string> algorithms{
			"labelGrid",
			"particlesToGrid",
			"extrapolateDepth2",
			"saveVelocityGrids",
			"applyBodyForces",
			"pressureSolve",
			"applyPressure",
			"gridToParticles",
			"extrapolateDepthN",
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

	delete m_particles;
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidSolver3D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	m_label = Mat2Di(m_gridWidth, m_gridHeight);
	m_p = Mat2Df(m_gridWidth, m_gridHeight);
	m_u = Mat2Df(m_gridWidth + 1, m_gridHeight);
	m_uSaved = Mat2Df(m_gridWidth + 1, m_gridHeight);
	m_v = Mat2Df(m_gridWidth, m_gridHeight + 1);
	m_vSaved = Mat2Df(m_gridWidth, m_gridHeight + 1);

	// init vel grids with unknown label value
	m_u.initValues(VEL_UNKNOWN);
	m_v.initValues(VEL_UNKNOWN);

	// read in initial geometry to populate label grid
	readInGeom2D(m_gridWidth, m_gridHeight, initialGeometryFile, m_label);

	// seed particles using label grid
	seedParticles(PARTICLES_PER_CELL, m_particles);
}

void FluidSolver3D::step() {
	// update the grid labels
	labelGrid();
	// transfer particle vel to grid
	particlesToGrid();
	// extrapolate fluid data out one cell for accurate divergence calculations
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, 2);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, 2);
	// save copy of current grid velocities for FLIP update
	saveVelocityGrids();
	// apply body forces on grid (gravity)
	applyBodyForces();
	// solve for pressure
	classicSolver solver(m_gridWidth, m_gridHeight, m_dx, m_dt, &m_label, &m_p, &m_u, &m_v);
	solver.pressureSolve();

	// apply pressure force
	applyPressure();
	// transfer grid velocities back to particles
	gridToParticles(PIC_WEIGHT);
	// advect particles
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridWidth);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridHeight);
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
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, 2);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, 2);
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
	classicSolver solver(m_gridWidth, m_gridHeight, m_dx, m_dt, &m_label, &m_p, &m_u, &m_v);
	solver.pressureSolve();
	m_timer->stop();

	// apply pressure force
	m_timer->start();
	applyPressure();
	m_timer->stop();
	// transfer grid velocities back to particles
	m_timer->start();
	gridToParticles(PIC_WEIGHT);
	m_timer->stop();
	// advect particles
	m_timer->start();
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridWidth);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridHeight);
	m_timer->stop();

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

void FluidSolver3D::saveLineData(std::ofstream *linesOut, float treshold) {
	if (linesOut->is_open()) {
		std::vector<glm::vec2> vertices = marchingSquares(m_p, m_gridWidth, m_gridHeight, treshold);
		for (int i = 0; i < vertices.size() - 1; i++) {
			(*linesOut) << vertices[i].x << " " << vertices[i].y << " ";
		}
		(*linesOut) << vertices[vertices.size() - 1].x << " " << vertices[vertices.size() - 1].y << "\n";
	}
}

void FluidSolver3D::saveLineData(std::ofstream *linesOut) {
	saveLineData(linesOut, SURFACE_THRESHOLD);
}

void FluidSolver3D::saveTriangleData(std::ofstream *vertOut, std::ofstream *indOut, std::ofstream *opacityOut, float treshold) {
	MarchingTrianglesData data = marchingSquaresTriangles(m_p, m_gridWidth, m_gridHeight, treshold);
	if (vertOut->is_open()) {
		for (int i = 0; i < data.vertices.size() - 1; i++) {
			(*vertOut) << data.vertices[i].x << " " << data.vertices[i].y << " ";
		}
		(*vertOut) << data.vertices[data.vertices.size() - 1].x << " " << data.vertices[data.vertices.size() - 1].y << "\n";
	}
	if (indOut->is_open()) {
		for (int i = 0; i < data.indices.size() - 1; i++) {
			(*indOut) << data.indices[i] << " ";
		}
		(*indOut) << data.indices[data.indices.size() - 1] << "\n";
	}
	if (opacityOut->is_open()) {
		for (int i = 0; i < data.opacities.size() - 1; i++) {
			(*opacityOut) << data.opacities[i] << " ";
		}
		(*opacityOut) << data.opacities[data.opacities.size() - 1] << "\n";
	}
}

void FluidSolver3D::saveTriangleData(std::ofstream *vertOut, std::ofstream *indOut, std::ofstream *opacityOut) {
	saveTriangleData(vertOut, indOut, opacityOut, SURFACE_THRESHOLD);
}

void FluidSolver3D::saveTimingData(std::ofstream *timingOut) {
	m_timer->writeTiming(timingOut);
}

std::vector<glm::vec2> FluidSolver3D::particleData() {
	std::vector<glm::vec2> particles;
	particles.reserve(m_particles->size());
	for (int i = 0; i < m_particles->size(); i++) {
		float x = 2 * m_particles->at(i).pos.x / (m_gridWidth * m_dx) - 1;
		float y = 2 * m_particles->at(i).pos.y / (m_gridHeight * m_dx) - 1;
		glm::vec2 p_i{ x,y };
		particles.push_back(p_i);
	}
	return particles;
}

std::vector<glm::vec2> FluidSolver3D::marchingSquares() {
	return marchingSquares(m_p, m_gridWidth, m_gridHeight, SURFACE_THRESHOLD);
}

MarchingTrianglesData FluidSolver3D::marchingSquaresTriangles() {
	return marchingSquaresTriangles(m_p, m_gridWidth, m_gridHeight, SURFACE_THRESHOLD);
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
void FluidSolver3D::seedParticles(int particlesPerCell, std::vector<Particle2D> *particleList) {
	// set random seed
	srand(time(NULL));
	// go through all cells marked fluid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label.get(i,j) == SimUtil::FLUID) {
				// seed randomly in 2x2 subgrid of the cell
				Vec2 cellCenter = getGridCellPosition(i, j, m_dx);
				Vec2 subCenters[] = {
					Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top left
					Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top right
					Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y - 0.25f*m_dx), // bottom right
					Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y - 0.25f*m_dx) // bottom left
				};
				// cycle through subgrid to place all particles
				for (int k = 0; k < particlesPerCell; k++) {
					// randomly jitter from subgrid center
					// give a random factor from [-0.24, 0.24] multiplied by dx
					float jitterX = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
					//jitterX = 0.0f; //DEBUG ALEX
					float jitterY = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
					//jitterY = 0.0f; //DEBUG ALEX
					Vec2 pos(subCenters[k % 4].x + jitterX, subCenters[k % 4].y + jitterY);
					Vec2 vel(0.0f, 0.0f);
					particleList->push_back(Particle2D(pos, vel));
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
			if (m_label.get(i,j) != SOLID) {
				m_label.set(i, j, AIR);
			}
		}
	}

	// mark any cell containing a particle FLUID
	for (int i = 0; i < m_particles->size(); i++) {
		// get cell containing the particle
		int *cell = getGridCellIndex(m_particles->at(i).pos, m_dx);
		m_label.set(cell[0], cell[1], FLUID);
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

	// For each component of velocity in each fluid grid cell
	// we calculate weighted average of particles around it defined
	// by a kernel function and set this as the vel in the cell. 

	// structures to accumulate grid numerator and denominator of weighted average before divide
	Mat2Dd uNum{ m_gridWidth + 1, m_gridHeight };
	Mat2Dd uDen{ m_gridWidth + 1, m_gridHeight };
	Mat2Dd vNum{ m_gridWidth, m_gridHeight + 1 };
	Mat2Dd vDen{ m_gridWidth, m_gridHeight + 1 };

	// clear accumulators
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			if (j < m_gridHeight) {
				uNum.set(i, j, 0.0);
				uDen.set(i, j, 0.0);
			}
			if (i < m_gridWidth) {
				vNum.set(i, j, 0.0);
				vDen.set(i, j, 0.0);
			}
		}
	}


	// loop over particles and accumulate num and den at each grid point
	if (USE_BSPLINE) {
		for (int p = 0; p < m_particles->size(); p++) {
			Particle2D curParticle = m_particles->at(p);
			int* ind = getGridCellIndex(curParticle.pos, m_dx);
			int i0 = ind[0];
			int j0 = ind[1];
			for (int i = i0 - 2; i <= i0 + 2; i++) {
				for (int j = j0 - 2; j <= j0 + 2; j++) {
					if (j < m_gridHeight) {
						double kernel = quadBSplineKernel(sub(curParticle.pos, getGridCellPosition(i - 0.5f, j, m_dx)));
						uNum.set(i, j, (uNum.get(i, j) + curParticle.vel.x * kernel));
						uDen.set(i, j, uDen.get(i, j) + kernel);
					}
					if (i < m_gridWidth) {
						double kernel = quadBSplineKernel(sub(curParticle.pos, getGridCellPosition(i, j - 0.5f, m_dx)));
						vNum.set(i, j, (vNum.get(i, j) + curParticle.vel.y * kernel));
						vDen.set(i, j, vDen.get(i, j) + kernel);
					}
				}
			}
		}
	}	
	else {
		for (int p = 0; p < m_particles->size(); p++) {
			Particle2D curParticle = m_particles->at(p);
			int* ind = getGridCellIndex(curParticle.pos, m_dx);
			int i0 = ind[0];
			int j0 = ind[1];
			for (int i = i0 - 1; i <= i0 + 1; i++) {
				for (int j = j0 - 1; j <= j0 + 1; j++) {
					if (j < m_gridHeight) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i - 0.5f, j, m_dx)));
						uNum.set(i, j, (uNum.get(i, j) + curParticle.vel.x * kernel));
						uDen.set(i, j, uDen.get(i, j) + kernel);
					}
					if (i < m_gridWidth) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i, j - 0.5f, m_dx)));
						vNum.set(i, j, (vNum.get(i, j) + curParticle.vel.y * kernel));
						vDen.set(i, j, vDen.get(i, j) + kernel);
					}
				}
			}
		}
	}

	// additional pass over grid to divide and update actual velocities
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			if (j < m_gridHeight) {
				if (uDen.get(i,j) != 0.0) {
					double val = uNum.get(i, j) / uDen.get(i, j);
					m_u.set(i, j, val);
				}
			}
			if (i < m_gridWidth) {
				if (vDen.get(i, j) != 0.0) {
					double val = vNum.get(i, j) / vDen.get(i, j);
					m_v.set(i, j, val);
				}
			}
		}
	}

	uNum.deleteGrid();
	uDen.deleteGrid();
	vNum.deleteGrid();
	vDen.deleteGrid();
}

/*
Extrapolates the data in fluid cells of the given grid out using a breadth-first
search technique.
Args:
grid - the grid with data to extrapolate
x, y - the grid dimensions
depth - the number of cells away from fluid cells to extrapolate to.
*/
void FluidSolver3D::extrapolateGridFluidData(Mat2Df &grid, int x, int y, int depth) {
	// initialize marker array 
	Mat2Di d{ x,y };
	// set d to 0 for known values, max int for unknown
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			if (grid.get(i,j) != VEL_UNKNOWN) {
				d.set(i, j, 0);
			} else {
				d.set(i, j, INT_MAX);
			}
		}
	}
	// define neighbors
	int numNeighbors = 8;
	int neighbors[8][2] = {
		{-1, 1}, // top left
		{-1, 0}, // middle left
		{-1, -1}, // bottom left
		{0, 1}, // top middle
		{0, -1}, // bottom middle
		{1, 1}, // top right
		{1, 0}, // middle right
		{1, -1} // bottom right
	};
	// initialize first wavefront
	std::vector<Vec2> W;
	int dim[2] = { x, y };
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			// current value is not known
			if (d.get(i,j) != 0) {
				int ind[2] = { i, j };
				if(hasNeighbors(d, dim, ind, neighbors, numNeighbors, 0)){
					// a neighbor is known
					d.set(i, j, 1);
					W.push_back(Vec2(i, j));
				}
			}
		}
	}
	// list of all wavefronts, only want to go through the given depth
	std::vector<std::vector<Vec2>> wavefronts;
	wavefronts.push_back(W);
	int curWave = 0;
	while (curWave < depth) {
		// get wavefront
		std::vector<Vec2> curW = wavefronts.at(curWave);
		// initialize next wavefront
		std::vector<Vec2> nextW;
		// go through current wave and extrapolate values
		for (int i = 0; i < curW.size(); i++) {
			Vec2 ind = curW.at(i);
			// average neighbors
			float avg = 0.0f;
			int numUsed = 0;
			for (int i = 0; i < numNeighbors; i++) {
				int offsetX = neighbors[i][0];
				int offsetY = neighbors[i][1];
				int neighborX = ind.x + offsetX;
				int neighborY = ind.y + offsetY;

				// make sure valid indices
				if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
					// only want to add to average if neighbor d is less than current d
					if (d.get(neighborX,neighborY) < d.get((int)ind.x,(int)ind.y)) {
						avg += grid.get(neighborX, neighborY);
						numUsed++;
					} else if (d.get(neighborX, neighborY) == INT_MAX) {
						d.set(neighborX, neighborY, d.get((int)ind.x, (int)ind.y) + 1);
						nextW.push_back(Vec2(neighborX, neighborY));
					}
				}
			}

			avg /= numUsed;
			// set current value to average of neighbors
			grid.set((int)ind.x, (int)ind.y, avg);
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
			m_uSaved.set(i, j, m_u.get(i, j));
		}
	}

	// save v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			m_vSaved.set(i, j, m_v.get(i, j));
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
			if (j < m_gridHeight) {
				// make sure we know the velocity
				if (m_u.get(i,j) != VEL_UNKNOWN) {
					// update u component
					m_u.set(i, j, m_u.get(i, j) + m_dt*GRAVITY.x);
				}
			}
			if (i < m_gridWidth) {
				if (m_v.get(i,j) != VEL_UNKNOWN) {
					// update v component
					m_v.set(i, j, m_v.get(i, j) + m_dt*GRAVITY.y);
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
			// update u
			if (i - 1 >= 0) {
				if (m_label.get(i - 1, j) == FLUID || m_label.get(i, j) == FLUID) {
					if (m_label.get(i - 1, j) == SOLID || m_label.get(i, j) == SOLID) {
						m_u.set(i, j, 0.0f); // usolid[i][j]
					} else {
						m_u.set(i, j, m_u.get(i, j) - scale * (m_p.get(i, j) - m_p.get(i - 1, j)));
					}
				} else {
					m_u.set(i, j, VEL_UNKNOWN);
				}
			} else {
				// edge of grid, keep the same velocity
			}

			// update v
			if (j - 1 >= 0) {
				if (m_label.get(i, j - 1) == FLUID || m_label.get(i, j) == FLUID) {
					if (m_label.get(i, j - 1) == SOLID || m_label.get(i, j) == SOLID) {
						m_v.set(i, j, 0.0f); // vsolid[i][j]
					}
					else {
						m_v.set(i, j, m_v.get(i, j) - scale * (m_p.get(i, j) - m_p.get(i, j - 1)));
					}
				} else {
					m_v.set(i, j, VEL_UNKNOWN);
				}
			} else {
				// edge of grid, keep the same velocity
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
	Mat2Df duGrid{ m_gridWidth + 1, m_gridHeight };
	Mat2Df dvGrid{ m_gridWidth, m_gridHeight + 1 };
	// calc u grid
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			duGrid.set(i, j, m_u.get(i, j) - m_uSaved.get(i, j));
		}
	}
	// calc v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			dvGrid.set(i, j, m_v.get(i, j) - m_vSaved.get(i, j));
		}
	}

	// go through particles and interpolate each velocity component
	// the update is a PIC/FLIP mix weighted with alpha
	// alpha = 1.0 is entirely PIC, alpha = 0.0 is all FLIP
	for (int i = 0; i < m_particles->size(); i++) {
		Particle2D *curParticle = &(m_particles->at(i));
		Vec2 picInterp = interpVel(m_u, m_v, curParticle->pos);
		Vec2 flipInterp = interpVel(duGrid, dvGrid, curParticle->pos);
		// u_new = alpha * interp(u_gridNew, x_p) + (1 - alpha) * (u_pOld + interp(u_dGrid, x_p))
		curParticle->vel = add(scale(picInterp, alpha), scale(add(curParticle->vel, flipInterp), 1.0f - alpha));
	}

	duGrid.deleteGrid();
	dvGrid.deleteGrid();
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
		Particle2D *curParticle = &(m_particles->at(i));
		float subTime = 0;
		bool finished = false;
		//float dT = m_dt / 4.999f;
		while (!finished) {
			Vec2 curVel = interpVel(m_u, m_v, curParticle->pos);

			// calc max substep size
			float dT = (C * m_dx) / (norm(curVel) + FLT_MIN);
			// update substep time so we don't go past normal time step
			if (subTime + dT >= m_dt) {
				dT = m_dt - subTime;
				finished = true;
			} else if (subTime + 2 * dT >= m_dt) {
				dT = 0.5f * (m_dt - subTime);
			}

			RK3(curParticle, curVel, dT, m_u, m_v);
			subTime += dT;

			if (curParticle->pos.x < 0 || curParticle->pos.y < 0 || isnan(curParticle->pos.x) || isnan(curParticle->pos.y)) {
				// there's been an error in RK3, just skip it
				std::cout << "RK3 error...skipping particle" << std::endl;
				break;
			}

			int *cell = getGridCellIndex(curParticle->pos, m_dx);
			int j = cell[0];
			int k = cell[1];
			if (m_label.get(j, k) == SOLID) {
				//std::cout << "Advected into SOLID, projecting back!\n";
				if (!projectParticle(curParticle, m_dx / 4.0f)) {
					std::cout << "RK3 error...skipping particle" << std::endl;
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
		int ind[2] = { cell[0], cell[1] };
		// if either of cells are negative or greater than sim dimensions it has left sim area
		if (ind[0] < 0 || ind[1] < 0 || ind[0] >= m_gridWidth || ind[1] >= m_gridHeight || isnan(m_particles->at(i).pos.x) || isnan(m_particles->at(i).pos.y)) {
			m_particles->erase(m_particles->begin() + i);
			numDeleted++;
			if (i >= m_particles->size()) {
				finished = true;
			}
		} else if (m_label.get(ind[0],ind[1]) == SOLID) {
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
std::vector<int> FluidSolver3D::checkNeighbors(Mat2Di &grid, int dim[2], int index[2], int neighbors[][2], int numNeighbors, int value) {
	std::vector<int> neighborsTrue;
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;

		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
			if (grid.get(neighborX, neighborY) == value) {
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
bool FluidSolver3D::hasNeighbors(Mat2Di &grid, int dim[2], int index[2], int neighbors[][2], int numNeighbors, int value) {
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;
		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
			if (grid.get(neighborX, neighborY) == value) {
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
double FluidSolver3D::trilinearHatKernel(SimUtil::Vec2 dist) {
	return hatFunction(dist.x / m_dx) * hatFunction(dist.y / m_dx);
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
Returns the value of the quadratic B-spline function for the
given distance (x, y).
*/
double FluidSolver3D::quadBSplineKernel(SimUtil::Vec2 dist) {
	return bSplineFunction(dist.x / m_dx) * bSplineFunction(dist.y / m_dx);
}

/*
Calculates the value of the bSpline function for given r
*/
double FluidSolver3D::bSplineFunction(double r) {
	double rAbs = abs(r);
	if (rAbs < 0.5) {
		return 0.75 - rAbs * rAbs;
	}
	if (rAbs < 1.5) {
		return 0.5 * (1.5 - rAbs)*(1.5 - rAbs);
	}
	return 0;
}

/*
Interpolates the value in the given velocity grid at the given position using bilinear interpolation.
Returns velocity unkown if position is not on simulation grid. 
Args:
uGrid - the u component grid to interpolate from
vGrid - the v component grid to interpolate from
pos - the position to interpolate at
*/
Vec2 FluidSolver3D::interpVel(SimUtil::Mat2Df &uGrid, SimUtil::Mat2Df &vGrid, Vec2 pos) {
	// get grid cell containing position
	int *cell = getGridCellIndex(pos, m_dx);
	int i = cell[0];
	int j = cell[1];
	// make sure this is a valid index
	if (i >= 0 && i < m_gridWidth && j >= 0 && j < m_gridHeight) {
		// get positions of u and v component stored on each side of cell
		Vec2 cellLoc = getGridCellPosition(i, j, m_dx);
		float offset = m_dx / 2.0f;
		float x1 = cellLoc.x - offset;
		float x2 = cellLoc.x + offset;
		float y1 = cellLoc.y - offset;
		float y2 = cellLoc.y + offset;
		// get actual values at these positions
		float u1 = uGrid.get(i, j);
		float u2 = uGrid.get(i + 1, j);
		float v1 = vGrid.get(i, j);
		float v2 = vGrid.get(i, j + 1);

		// the interpolated values
		float u = ((x2 - pos.x) / (x2 - x1)) * u1 + ((pos.x - x1) / (x2 - x1)) * u2;
		float v = ((y2 - pos.y) / (y2 - y1)) * v1 + ((pos.y - y1) / (y2 - y1)) * v2;
		return Vec2(u, v);
	} else {
		return Vec2(VEL_UNKNOWN, VEL_UNKNOWN);
	}
}

/*
Advects a particle using Runge-Kutta 3 through the given velocity field.
Args:
particle - the particle to advect
initVel - the particles initial velocity in the current field, can leave UNKNOWN
dt - the time step
uGrid/vGrid - the velocity grids to advect through
*/
void FluidSolver3D::RK3(SimUtil::Particle2D *particle, SimUtil::Vec2 initVel, float dt, SimUtil::Mat2Df &uGrid, SimUtil::Mat2Df &vGrid) {
	if (initVel.x == VEL_UNKNOWN && initVel.y == VEL_UNKNOWN) {
		initVel = interpVel(uGrid, vGrid, particle->pos);
	}

	Vec2 k1 = initVel;
	Vec2 k2 = interpVel(uGrid, vGrid, add(particle->pos, scale(k1, 0.5f*dt)));
	Vec2 k3 = interpVel(uGrid, vGrid, add(particle->pos, scale(k2, 0.75f*dt)));
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

bool FluidSolver3D::projectParticle(Particle2D *particle, float dx) {
	// project back into fluid
	// find neighbors that are fluid
	// define neighbors
	int numNeighbors = 8;
	int neighbors[8][2] = {
		{ -1, 1 }, // top left
		{ -1, 0 }, // middle left
		{ -1, -1 }, // bottom left
		{ 0, 1 }, // top middle
		{ 0, -1 }, // bottom middle
		{ 1, 1 }, // top right
		{ 1, 0 }, // middle right
		{ 1, -1 } // bottom right
	};
	int dim[2] = { m_gridWidth, m_gridHeight };
	int *cell = getGridCellIndex(particle->pos, m_dx);
	int index[2] = { cell[0], cell[1] };
	// get neighbors that are fluid
	std::vector<int> neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, FLUID);
	if (neighborInd.size() == 0) {
		// try with air
		neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, AIR);
	}
	// find closest to particle
	int closestInd = -1;
	float closestDist = std::numeric_limits<float>::max();
	Vec2 closestVec(0.0f, 0.0f);
	for (int j = 0; j < neighborInd.size(); j++) {
		// get vec from particle to neighbor ind
		int ind[2] = { cell[0] + neighbors[neighborInd.at(j)][0], cell[1] + neighbors[neighborInd.at(j)][1] };
		Vec2 cellPos = getGridCellPosition(ind[0], ind[1], m_dx);
		Vec2 distVec = sub(cellPos, particle->pos);
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
		Vec2 projectVec(0.0f, 0.0f);
		if (closestInd == 1) { // middle left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 3) { // top middle
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 4) { // bottom middle
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 6) { // middle right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 5) { // top right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 0) { // top left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 2) { // bottom left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 7) { // bottom right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}

		particle->pos = add(particle->pos, projectVec);

		return true;
	}
}

/*
Returns the line data of the isocontour in a given grid as a vertices.
each line has a startpoint and an end point stored as a vec2 in the vector
Args:
grid - given grid to extract the isocontour
width - gridwidth
height - gridheigt
tol - threshold of for the isocontour
*/
std::vector<glm::vec2> FluidSolver3D::marchingSquares(SimUtil::Mat2Df &grid, int width, int height, float tol) {
	std::vector<glm::vec2> points;

	//16 different cases for one cell for the marching squares
	std::vector<std::vector<float>> cases;
	cases.push_back(std::vector<float> {});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 1.0f});
	cases.push_back(std::vector<float> {0.5f, 1.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.5f, 0.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 1.0f, 0.5f, 0.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.5f, 0.0f, 0.5f, 1.0f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 0.0f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 0.0f});
	cases.push_back(std::vector<float> {0.5f, 1.0f, 0.5f, 0.0f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 0.0f, 0.5f, 1.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.5f, 0.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.5f, 1.0f, 1.0f, 0.5f});
	cases.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 1.0f});
	cases.push_back(std::vector<float> {});

	for (int i = 0; i < width - 1; i++) {
		for (int j = 0; j < height - 1; j++) {
			//determine which of the 16 different cases for one square exists
			int selectCase = 0;
			if (grid.get(i, j) > tol)
				selectCase += 8;
			if (grid.get(i + 1, j) > tol)
				selectCase += 4;
			if (grid.get(i + 1, j + 1) > tol)
				selectCase += 2;
			if (grid.get(i, j + 1) > tol)
				selectCase += 1;
			for (int k = 0; k < cases[selectCase].size() / 2; k++) {
				float offsetX = cases[selectCase][2 * k];
				float offsetY = cases[selectCase][2 * k + 1];

				//the offset of 0.5f indicates that this axis needs interpolation
				if (offsetX == 0.5f) {
					int offsetJ = (int)offsetY;
					offsetX = 1 / (grid.get(i + 1, j + offsetJ) - grid.get(i, j + offsetJ)) * (tol - grid.get(i, j + offsetJ));
				}
				else {
					int offsetI = (int)offsetX;
					offsetY = 1 / (grid.get(i + offsetI, j + 1) - grid.get(i + offsetI, j)) * (tol - grid.get(i + offsetI, j));
				}

				double x = 2 * (i + offsetX) / (width - 1) - 1;
				double y = 2 * (j + offsetY) / (height - 1) - 1;

				points.push_back(glm::vec2{ x,y });
			}
		}
	}
	return points;
}

/*
Returns the triangle data of the isocontour in a given grid as a struct of 3 vectors.
first containing vertices as vec2
second containing vertex indices for index buffering
third containing opacity values per vertex for pressure dependend opacity
Args:
grid - given grid to extract the isocontour
width - gridwidth
height - gridheigt
tol - threshold of for the isocontour
*/
MarchingTrianglesData FluidSolver3D::marchingSquaresTriangles(SimUtil::Mat2Df &grid, int width, int height, float tol) {
	std::vector<glm::vec2> vertices;
	std::vector<int> globalIndices;
	std::vector<float> opacities;
	//stores the maximum pressure for interpolation
	float pMax = grid.max();
	int curInd = 0;
	//16 different cases for one cell for the marching squares
	std::vector<std::vector<float>> points;
	std::vector<std::vector<int>> indices;
	points.push_back(std::vector<float> {});
	indices.push_back(std::vector<int> {});
	points.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 1.0f, 0.0f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2});
	points.push_back(std::vector<float> {0.5f, 1.0f, 1.0f, 0.5f, 1.0f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2});
	points.push_back(std::vector<float> {0.0f, 0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.5f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 0, 3});
	points.push_back(std::vector<float> {0.5f, 0.0f, 1.0f, 0.5f, 1.0f, 0.0f});
	indices.push_back(std::vector<int> {0, 1, 2});
	points.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 1.0f, 0.0f, 1.0f, 0.5f, 0.0f, 1.0f, 0.5f, 1.0f, 0.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 3, 4, 5});
	points.push_back(std::vector<float> {0.5f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.5f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 3, 0});
	points.push_back(std::vector<float> {0.0f, 0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.5f, 0.0f, 1.0f, 0.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 0, 2, 3, 3, 2, 4});
	points.push_back(std::vector<float> {0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f});
	indices.push_back(std::vector<int> {0, 1, 2});
	points.push_back(std::vector<float> {0.0f, 0.0f, 0.5f, 0.0f, 0.0f, 1.0f, 0.5f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 3, 1});
	points.push_back(std::vector<float> {0.0f, 0.0f, 0.5f, 0.0f, 0.0f, 0.5f, 1.0f, 1.0f, 1.0f, 0.5f, 0.5f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 3, 4, 5});
	points.push_back(std::vector<float> {0.0f, 0.0f, 0.5f, 0.0f, 0.0f, 1.0f, 1.0f, 0.5f, 1.0f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 1, 2, 3, 3, 2, 4});
	points.push_back(std::vector<float> {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 1.0f, 0.5f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 3, 1});
	points.push_back(std::vector<float> {0.0f, 0.0f, 0.0f, 1.0f, 0.5f, 1.0f, 1.0f, 0.5f, 1.0f, 0.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 0, 2, 3, 0, 4, 3});
	points.push_back(std::vector<float> {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 1, 3, 3, 1, 4});
	points.push_back(std::vector<float> {0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f});
	indices.push_back(std::vector<int> {0, 1, 2, 2, 1, 3});

	for (int i = 0; i < width - 1; i++) {
		for (int j = 0; j < height - 1; j++) {
			//determine which of the 16 different cases for one square exists
			int selectCase = 0;
			if (grid.get(i, j) > tol)
				selectCase += 8;
			if (grid.get(i + 1, j) > tol)
				selectCase += 4;
			if (grid.get(i + 1, j + 1) > tol)
				selectCase += 2;
			if (grid.get(i, j + 1) > tol)
				selectCase += 1;
			for (int k = 0; k < indices[selectCase].size(); k++) {
				globalIndices.push_back(indices[selectCase][k] + curInd);
			}
			for (int k = 0; k < points[selectCase].size() / 2; k++) {
				float offsetX = points[selectCase][2 * k];
				float offsetY = points[selectCase][2 * k + 1];

				float p = grid.get((int)(i + offsetX), (int)(j + offsetY));

				//the offset of 0.5f indicates that this axis needs interpolation
				if (offsetX == 0.5f) {
					int offsetJ = (int)offsetY;
					offsetX = 1 / (grid.get(i + 1, j + offsetJ) - grid.get(i, j + offsetJ)) * (tol - grid.get(i, j + offsetJ));

					p = grid.get(i, j + offsetJ) + (grid.get(i + 1, j + offsetJ) - grid.get(i, j + offsetJ)) * offsetX;
				}
				if (offsetY == 0.5f) {
					int offsetI = (int)offsetX;
					offsetY = 1 / (grid.get(i + offsetI, j + 1) - grid.get(i + offsetI, j)) * (tol - grid.get(i + offsetI, j));

					p = grid.get(i + offsetI, j) + (grid.get(i + offsetI, j + 1) - grid.get(i + offsetI, j)) * offsetY;
				}

				float x = 2 * (i + offsetX) / (width - 1) - 1;
				float y = 2 * (j + offsetY) / (height - 1) - 1;

				//1 for maximum value in grid, 0 for tolerance
				float opacity = (p - tol) / (pMax - tol);

				vertices.push_back(glm::vec2{ x,y });
				opacities.push_back(opacity);
				//increment the current Indice value
				curInd++;
			}

		}
	}
	return MarchingTrianglesData(vertices, globalIndices, opacities);
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

