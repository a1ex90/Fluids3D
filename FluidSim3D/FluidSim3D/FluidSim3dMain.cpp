//----------------------------------------------------------------------
// FluidSolver Project
// by Alex Sommer
//----------------------------------------------------------------------

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>

#include "FluidSolver3d.h"
#include "FluidRenderer3d.h"
#include "SimUtil.h"
#include "MarchingCubes.h"

//----------------------------------------------------------------------
// Execution Options
//----------------------------------------------------------------------

// whether to run the simulation
const bool RUN_SIM = false;
// whether to run rendering
const bool RUN_RENDERING = true;
// wheater to do a realtime simulation
const bool REALTIME_SIM = false;
/* 
choose a visualization mode for rendering
0 - Draw Surface as Line
1 - Draw Fluid as Triangle Mesh
2 - Draw Fluid as pressure dependend opacity Triangle Mesh
3 - Draw Fluid as Particles
4 - Draw Fluid as Particles + Lines
*/
const int VISUALIZATION_MODE = 3;

//----------------------------------------------------------------------
// Simulation Parameters
//----------------------------------------------------------------------

// resolution of the grid to use (width, height)
const int GRID_WIDTH = 25;
const int GRID_HEIGHT = 25;
const int GRID_DEPTH = 25;
// grid cell width (in meters)
const float GRID_CELL_WIDTH = 0.005f;
// simulation time step (in seconds)
const float TIME_STEP = 0.01f;

//----------------------------------------------------------------------
// I/O Parameters
//----------------------------------------------------------------------

// input file for initial system state - grid marked solid, fluid, or air
const std::string INITIAL_GEOMETRY_FILE_IN = "geo_small3.txt";
// output file for watersurface line data
const std::string DATA_FILES_OUT = "shit";


//----------------------------------------------------------------------
// Global Variables
//----------------------------------------------------------------------

// the number of frames to simulate
const int NUM_SIM_FRAMES = 100;
// frame rate for render (fps)
const float FRAME_RATE = 25.0f;
// time step between outputted frames
const float FRAME_TIME_STEP = 1.0f / FRAME_RATE;

//----------------------------------------------------------------------
// Main Function
//----------------------------------------------------------------------

int main(int argc, char** argv) {
	if (REALTIME_SIM) {
		FluidSolver3D solver(GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH, GRID_CELL_WIDTH, TIME_STEP);
		solver.init(INITIAL_GEOMETRY_FILE_IN);
		//FluidRenderer3D render{ INITIAL_GEOMETRY_FILE_IN, VISUALIZATION_MODE };
		FluidRenderer3D render( solver.getGeometry(),GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH, VISUALIZATION_MODE );
		for (int frame = 0; frame < 4*NUM_SIM_FRAMES; frame++) {
			auto start = std::chrono::system_clock::now();
			solver.step();
			if (VISUALIZATION_MODE == 0) {
				//render.draw(solver.marchingSquares());
			}
			else if (VISUALIZATION_MODE == 3) {
				render.drawP(solver.particleData());
			}
			else {
				/*SimUtil::MarchingTrianglesData data = solver.marchingSquaresTriangles();
				if (VISUALIZATION_MODE == 1) {
					render.draw(data.vertices, data.indices);
				}
				else if (VISUALIZATION_MODE == 2) {
					render.draw(data.vertices, data.indices, data.opacities);
				}*/
			}
			if (frame == 0) {
				//Wait for keypress to start
				std::cin.ignore();
			}
			bool sleep = true;
			while (sleep) {
				auto now = std::chrono::system_clock::now();
				auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
				if (elapsed.count() > 1000*TIME_STEP) {
					sleep = false;
				}
			}
		}
	}

	/*if (RUN_SIM) {
		// open and clear output file
		std::ofstream *particlesOut = new std::ofstream(DATA_FILES_OUT + "-part.csv", std::ofstream::trunc);
		std::ofstream *linesOut = new std::ofstream(DATA_FILES_OUT + "-lines.csv", std::ofstream::trunc);
		std::ofstream *timingOut = new std::ofstream(DATA_FILES_OUT + "-timing.csv", std::ofstream::trunc);
		std::ofstream *trianglesVert = new std::ofstream(DATA_FILES_OUT + "-vert.csv", std::ofstream::trunc);
		std::ofstream *trianglesInd = new std::ofstream(DATA_FILES_OUT + "-ind.csv", std::ofstream::trunc);
		std::ofstream *trianglesOpa = new std::ofstream(DATA_FILES_OUT + "-opacity.csv", std::ofstream::trunc);

		FluidSolver3D solver(GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH, TIME_STEP);
		std::cout << "Simulating Frame 1" << std::endl;
		solver.init(INITIAL_GEOMETRY_FILE_IN);
		
		// run simulation
		// save initial frame
		//ADD INITIAL STEP HERE FOR MARCHING SQUARES TO WORK
		solver.step();
		solver.saveParticleData(particlesOut);
		solver.saveLineData(linesOut);
		solver.saveTriangleData(trianglesVert, trianglesInd, trianglesOpa);
		float t = 0.0f;
		for (int framesOut = 1; framesOut < NUM_SIM_FRAMES; framesOut++) {
			std::cout << "Simulating Frame " << framesOut + 1 << std::endl;
			
			float subTime = 0.0f;
			while (subTime < FRAME_TIME_STEP) {
				// perform sim time step
				//solver.step();
				solver.stepTiming();

				subTime += TIME_STEP;
				t += TIME_STEP;
			}
			// render at current time
			solver.saveParticleData(particlesOut);
			solver.saveLineData(linesOut);
			solver.saveTriangleData(trianglesVert, trianglesInd, trianglesOpa);
			
		}
		solver.saveTimingData(timingOut);
		// cleanup
		particlesOut->close();
		linesOut->close();
		trianglesVert->close();
		trianglesInd->close();
		trianglesOpa->close();
		timingOut->close();
		delete particlesOut;
		delete linesOut;
		delete trianglesVert;
		delete trianglesInd;
		delete trianglesOpa;
		delete timingOut;
	}*/

	if (RUN_RENDERING) {
		int CASENO = 83;
		SimUtil::Mat3Di empty{ 2,2,2 };
		empty.initValues(SimUtil::AIR);
		SimUtil::Mat3Df caseMat{ 2,2,2 };
		MarchingCubes::initCase(caseMat, CASENO);
		std::vector<std::vector<glm::vec3>> cubeCases; 
		std::vector<std::vector<int>> cubeIndices;
		
		std::vector<glm::vec3> darkDots;
		std::vector<glm::vec3> brightDots;
		MarchingCubes::corners(darkDots, brightDots, CASENO);
		MarchingCubes::initMarchingCubesCases(cubeCases, cubeIndices);
		int count = 0;
		for (int i = 0; i < 256; i++) {
			
			if (cubeCases[i].empty()) {
				count++;
			}
			
		}
		std::cout << count << std::endl;
		SimUtil::Mesh3D caseMesh = MarchingCubes::meshData(caseMat, cubeCases, cubeIndices, 2, 2, 2, 0.0f);


		FluidRenderer3D render(&empty, 2, 2, 2, VISUALIZATION_MODE);
		while (true) {
			render.drawCubes(caseMesh.vertices, caseMesh.indices, darkDots, brightDots);
		}
		
	}

	return 0;
}