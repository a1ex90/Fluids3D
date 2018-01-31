#include "SimUtil.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

namespace SimUtil {
	//----------------------------------------------------------------------
	// Mat3Di Class
	//----------------------------------------------------------------------
	Mat3Di::Mat3Di() {
		m_width = 0;
		m_height = 0;
		m_depth = 0;
		m_grid = 0;
	}
	Mat3Di::Mat3Di(int width, int height, int depth) {
		m_width = width;
		m_height = height;
		m_depth = depth;
		m_grid = new int[width * height * depth];
	}

	void Mat3Di::initValues(int val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					m_grid[i + m_width * j + k * m_width * m_height] = val;
				}		
			}
		}
	}

	Mat3Di::~Mat3Di() {
	}

	void Mat3Di::set(int i, int j, int k, int val) {
		m_grid[i + m_width * j + k * m_width * m_height] = val;
	}

	int Mat3Di::get(int i, int j, int k) {
		return m_grid[i + m_width * j + k * m_width * m_height];
	}

	int Mat3Di::height() {
		return m_height;
	}

	int Mat3Di::width() {
		return m_width;
	}

	int Mat3Di::depth() {
		return m_depth;
	}

	int Mat3Di::max() {
		int maxVal = std::numeric_limits<int>::lowest();
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					if (m_grid[i + m_width * j + k * m_width * m_height] > maxVal) {
						maxVal = m_grid[i + m_width * j + k * m_width * m_height];
					}
				}		
			}
		}
		return maxVal;
	}

	double Mat3Di::dot(Mat3Di o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					dotProd += m_grid[i + m_width * j + k * m_width * m_height] * o_mat.get(i, j, k);
				}
			}
		}
		return dotProd;
	}

	void Mat3Di::deleteGrid() {
		m_height = 0;
		m_width = 0;
		m_depth = 0;
		delete[] m_grid;
	}

	//----------------------------------------------------------------------
	// Mat3Df Class
	//----------------------------------------------------------------------
	Mat3Df::Mat3Df() {
		m_width = 0;
		m_height = 0;
		m_grid = 0;
	}
	Mat3Df::Mat3Df(int width, int height, int depth) {
		m_width = width;
		m_height = height;
		m_depth = depth;
		m_grid = new float[width * height * depth];
	}

	void Mat3Df::initValues(float val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					m_grid[i + m_width * j + k * m_width * m_height] = val;
				}
			}
		}
	}

	Mat3Df::~Mat3Df() {
	}

	void Mat3Df::set(int i, int j, int k, float val) {
		m_grid[i + m_width * j + k * m_width * m_height] = val;
	}

	float Mat3Df::get(int i, int j, int k) {
		return m_grid[i + m_width * j + k * m_width * m_height];
	}

	int Mat3Df::height() {
		return m_height;
	}

	int Mat3Df::width() {
		return m_width;
	}

	int Mat3Df::depth() {
		return m_depth;
	}

	float Mat3Df::max() {
		float maxVal = std::numeric_limits<float>::lowest();
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				if (m_grid[i + m_width * j] > maxVal) {
					maxVal = m_grid[i + m_width * j];
				}
			}
		}
		return maxVal;
	}

	double Mat3Df::dot(Mat3Df o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					dotProd += m_grid[i + m_width * j + k * m_width * m_height] * o_mat.get(i, j, k);
				}
			}
		}
		return dotProd;
	}

	void Mat3Df::deleteGrid() {
		m_height = 0;
		m_width = 0;
		m_depth = 0;
		delete[] m_grid;
	}

	//----------------------------------------------------------------------
	// Mat3Dd Class
	//----------------------------------------------------------------------
	Mat3Dd::Mat3Dd() {
		m_width = 0;
		m_height = 0;
		m_grid = 0;
	}
	Mat3Dd::Mat3Dd(int width, int height, int depth) {
		m_width = width;
		m_height = height;
		m_depth = depth;
		m_grid = new double[width * height * depth];
	}

	void Mat3Dd::initValues(double val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					m_grid[i + m_width * j + k * m_width * m_height] = val;
				}
			}
		}
	}

	Mat3Dd::~Mat3Dd() {
	}

	void Mat3Dd::set(int i, int j, int k, double val) {
		m_grid[i + m_width * j + k * m_width * m_height] = val;
	}

	double Mat3Dd::get(int i, int j, int k) {
		return m_grid[i + m_width * j + k * m_width * m_height];
	}

	int Mat3Dd::height() {
		return m_height;
	}

	int Mat3Dd::width() {
		return m_width;
	}

	int Mat3Dd::depth() {
		return m_depth;
	}

	double Mat3Dd::max() {
		double maxVal = std::numeric_limits<double>::lowest();
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				if (m_grid[i + m_width * j] > maxVal) {
					maxVal = m_grid[i + m_width * j];
				}
			}
		}
		return maxVal;
	}

	double Mat3Dd::dot(Mat3Dd o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				for (int k = 0; k < m_depth; k++) {
					dotProd += m_grid[i + m_width * j + k * m_width * m_height] * o_mat.get(i, j, k);
				}
			}
		}
		return dotProd;
	}

	double* Mat3Dd::data() {
		return m_grid;
	}

	void Mat3Dd::deleteGrid() {
		m_height = 0;
		m_width = 0;
		m_depth = 0;
		delete[] m_grid;
	}

	void readInGeom3D(int x, int y, int z, std::string geomFileName, Mat3Di &grid) {
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				grid.set(i, j, 0, SimUtil::SOLID);
				grid.set(i, j, 1, SimUtil::SOLID);
				grid.set(i, j, 2, SimUtil::SOLID);
				grid.set(i, j, z - 3, SimUtil::SOLID);
				grid.set(i, j, z - 2, SimUtil::SOLID);
				grid.set(i, j, z - 1, SimUtil::SOLID);
			}
		}
		for (int k = 3; k < z - 3; k++) {
			readInGeom2D(x, y, k, geomFileName, grid);
		}	
	}

	void readInGeom2D(int x, int y, int fixedZ, std::string geomFileName, Mat3Di &grid) {
		//TODO FIND SOME READ IN
		// open the geometry file
		std::ifstream geomFile(geomFileName);
		if (geomFile.is_open()) {
			std::string lineStr;
			// parse file based on given dimensions, will error if file does not match these
			// fills grid so that [0][0] is at bottom left corner of simulation
			for (int i = y - 1; i >= 0; i--) {
				std::getline(geomFile, lineStr);
				for (int j = 0; j < x; j++) {
					switch (lineStr[j]) {
					case 'f':
						grid.set(j, i, fixedZ, SimUtil::FLUID);
						break;
					case 's':
						grid.set(j, i, fixedZ, SimUtil::SOLID);
						break;
					case 'a':
						grid.set(j, i, fixedZ, SimUtil::AIR);
						break;
					}
				}
			}
			geomFile.close();
		}
	}

	Vec3 getGridCellPosition(float i, float j, float k, float dx) {
		return Vec3(i*dx + 0.5f*dx, j*dx + 0.5f*dx, k*dx + 0.0f*dx);
	}

	int* getGridCellIndex(Vec3 pos, float dx) {
		int index[3] = { (int)(pos.x / dx), (int)(pos.y / dx), (int)(pos.z / dx) };
		return index;
	}

	Vec3 add(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
	}

	Vec3 sub(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
	}

	Vec3 scale(Vec3 vec1, float scalar) {
		return Vec3(scalar * vec1.x, scalar * vec1.y, scalar * vec1.z);
	}

	float norm(Vec3 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}

}