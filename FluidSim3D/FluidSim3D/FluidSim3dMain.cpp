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

// whether to run rendering
const bool DISPLAY_MARCHINGCUBES = false;
// wheater to do a realtime simulation
const bool REALTIME_SIM = true;
// wheater to manipulate gravity by user input
const bool MANIPULATION = true;
/* 
choose a visualization mode for rendering
0 - Draw Fluid as Particles
1 - Draw Fluid as Triangle Mesh
*/
const int VISUALIZATION_MODE = 1;

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
const std::string INITIAL_GEOMETRY_FILE_IN = "geo_small4.txt";
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
		FluidRenderer3D render( solver.getGeometry(),GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH );
		/*for (int frame = 0; frame < 4*NUM_SIM_FRAMES; frame++) {
			auto start = std::chrono::system_clock::now();
			solver.step();
			if (VISUALIZATION_MODE == 0) {
				render.drawP(solver.particleData());
			}
			else if (VISUALIZATION_MODE == 1) {
				SimUtil::Mesh3D data = solver.meshData();
				render.draw(data.vertices, data.normals, data.indices);
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
		}*/
		solver.step();
		SimUtil::Mesh3D data = solver.meshData();
		auto start = std::chrono::system_clock::now();
		bool newFrame = true;
		while (!render.isClosed()) {
			render.draw(data.vertices, data.normals, data.indices);
			if (!render.isPaused() && newFrame) {
				start = std::chrono::system_clock::now();
				if (MANIPULATION) {
					solver.updateOrientation(render.currentOrientation());
				}
				solver.step();
				data = solver.meshData();
				newFrame = false;
			}
			if (render.isPaused() && render.forwardPressed()) {
				if (MANIPULATION) {
					solver.updateOrientation(render.currentOrientation());
				}
				solver.step();
				data = solver.meshData();
			}
			//Check if it's time for a new frame
			if (!render.isPaused() && !newFrame) {
				auto now = std::chrono::system_clock::now();
				auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
				if (elapsed.count() > 1000 * TIME_STEP) {
					newFrame = true;
				}
			}
		}
	}
	if (DISPLAY_MARCHINGCUBES) {
		int CASENO = 139;
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

		SimUtil::Mesh3D caseMesh = MarchingCubes::meshData(caseMat, cubeCases, cubeIndices, 2, 2, 2, 0.0f);


		FluidRenderer3D render(&empty, 2, 2, 2);
		while (true) {
			render.drawCubes(caseMesh.vertices, caseMesh.indices, darkDots, brightDots);
		}
		
	}

	return 0;
}