#include "SimUtil.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

namespace SimUtil {
	//----------------------------------------------------------------------
	// Mat2Di Class
	//----------------------------------------------------------------------
	Mat2Di::Mat2Di() {
		m_width = 0;
		m_height = 0;
		m_grid = 0;
	}
	Mat2Di::Mat2Di(int width, int height) {
		m_width = width;
		m_height = height;
		m_grid = new int[width * height];
	}

	void Mat2Di::initValues(int val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				m_grid[i + m_width * j] = val;
			}
		}
	}

	Mat2Di::~Mat2Di() {
	}

	void Mat2Di::set(int i, int j, int val) {
		m_grid[i + m_width * j] = val;
	}

	int Mat2Di::get(int i, int j) {
		return m_grid[i + m_width * j];
	}

	int Mat2Di::height() {
		return m_height;
	}

	int Mat2Di::width() {
		return m_width;
	}

	int Mat2Di::max() {
		int maxVal = std::numeric_limits<int>::lowest();
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				if (m_grid[i + m_width * j] > maxVal) {
					maxVal = m_grid[i + m_width * j];
				}
			}
		}
		return maxVal;
	}

	double Mat2Di::dot(Mat2Di o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				dotProd += m_grid[i + m_width * j] * o_mat.get(i,j);
			}
		}
		return dotProd;
	}

	void Mat2Di::deleteGrid() {
		m_height = 0;
		m_width = 0;
		delete[] m_grid;
	}

	//----------------------------------------------------------------------
	// Mat2Df Class
	//----------------------------------------------------------------------
	Mat2Df::Mat2Df() {
		m_width = 0;
		m_height = 0;
		m_grid = 0;
	}
	Mat2Df::Mat2Df(int width, int height) {
		m_width = width;
		m_height = height;
		m_grid = new float[width * height];
	}

	void Mat2Df::initValues(float val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				m_grid[i + m_width * j] = val;
			}
		}
	}

	Mat2Df::~Mat2Df() {
	}

	void Mat2Df::set(int i, int j, float val) {
		m_grid[i + m_width * j] = val;
	}

	float Mat2Df::get(int i, int j) {
		return m_grid[i + m_width * j];
	}

	int Mat2Df::height() {
		return m_height;
	}

	int Mat2Df::width() {
		return m_width;
	}

	float Mat2Df::max() {
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

	double Mat2Df::dot(Mat2Df o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				dotProd += m_grid[i + m_width * j] * o_mat.get(i, j);
			}
		}
		return dotProd;
	}

	void Mat2Df::deleteGrid() {
		m_height = 0;
		m_width = 0;
		delete[] m_grid;
	}

	//----------------------------------------------------------------------
	// Mat2Dd Class
	//----------------------------------------------------------------------
	Mat2Dd::Mat2Dd() {
		m_width = 0;
		m_height = 0;
		m_grid = 0;
	}
	Mat2Dd::Mat2Dd(int width, int height) {
		m_width = width;
		m_height = height;
		m_grid = new double[width * height];
	}

	void Mat2Dd::initValues(double val) {
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				m_grid[i + m_width * j] = val;
			}
		}
	}

	Mat2Dd::~Mat2Dd() {
	}

	void Mat2Dd::set(int i, int j, double val) {
		m_grid[i + m_width * j] = val;
	}

	double Mat2Dd::get(int i, int j) {
		return m_grid[i + m_width * j];
	}

	int Mat2Dd::height() {
		return m_height;
	}

	int Mat2Dd::width() {
		return m_width;
	}

	double Mat2Dd::max() {
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

	double Mat2Dd::dot(Mat2Dd o_mat) {
		double dotProd = 0.0;
		for (int i = 0; i < m_width; i++) {
			for (int j = 0; j < m_height; j++) {
				dotProd += m_grid[i + m_width * j] * o_mat.get(i, j);
			}
		}
		return dotProd;
	}

	double* Mat2Dd::data() {
		return m_grid;
	}

	void Mat2Dd::deleteGrid() {
		m_height = 0;
		m_width = 0;
		delete[] m_grid;
	}

	void readInGeom2D(int x, int y, std::string geomFileName, Mat2Di &grid) {
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
						grid.set(j, i, SimUtil::FLUID);
						break;
					case 's':
						grid.set(j, i, SimUtil::SOLID);
						break;
					case 'a':
						grid.set(j, i, SimUtil::AIR);
						break;
					}
				}
			}
			geomFile.close();
		}
	}

	Vec2 getGridCellPosition(float i, float j, float dx) {
		return Vec2(i*dx + 0.5f*dx, j*dx + 0.5f*dx);
	}

	int* getGridCellIndex(Vec2 pos, float dx) {
		int index[2] = { (int)(pos.x / dx), (int)(pos.y / dx) };
		return index;
	}

	Vec2 add(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x + vec2.x, vec1.y + vec2.y);
	}

	Vec2 sub(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x - vec2.x, vec1.y - vec2.y);
	}

	Vec2 scale(Vec2 vec1, float scalar) {
		return Vec2(scalar * vec1.x, scalar * vec1.y);
	}

	float norm(Vec2 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y);
	}

}