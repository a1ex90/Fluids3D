#include "SimUtil.h"

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

	void readInGeom3D(int x, int y, int z, int borderCount, int emptyFrontCount, int emptyBackCount, int zDepth, std::string geomFileName, Mat3Di &grid) {
		//init borders
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				for (int k = 0; k < borderCount; k++){
					grid.set(i, j, k, SimUtil::SOLID);
					grid.set(i, j, z - (k + 1), SimUtil::SOLID);
				}
			}
		}
		for (int k = borderCount; k < z - borderCount; k++) {
			for (int i = 0; i < x; i++) {
				for (int j = 0; j < borderCount; j++) {
					grid.set(i, j, k, SimUtil::SOLID);
					grid.set(i, y - (j + 1), k, SimUtil::SOLID);
				}
			}
			for (int i = 0; i < borderCount; i++) {
				for (int j = borderCount; j < y - borderCount; j++) {
					grid.set(i, j, k, SimUtil::SOLID);
					grid.set(x - (i + 1), j, k, SimUtil::SOLID);
				}
			}
		}
		//init empty layers
		for (int i = borderCount; i < x - borderCount; i++) {
			for (int j = borderCount; j < y - borderCount; j++) {
				for (int k = borderCount; k < borderCount + emptyFrontCount; k++) {
					grid.set(i, j, k, SimUtil::AIR);
				}
				for (int k = z - (borderCount + 1); k > z - (borderCount + emptyBackCount + 1); k--) {
					grid.set(i, j, k, SimUtil::AIR);
				}
			}
		}
		//fill the rest of the scene
		for (int k = borderCount + emptyFrontCount; k < z - (borderCount + emptyBackCount); k++) {
			readInGeom2D(x, y, k, borderCount, geomFileName, grid);
		}
	}

	void readInGeom2D(int x, int y, int fixedZ, int borderCount, std::string geomFileName, Mat3Di &grid) {
		//TODO FIND SOME READ IN
		// open the geometry file
		std::ifstream geomFile(geomFileName);
		if (geomFile.is_open()) {
			std::string lineStr;
			//Skip the first 4 files since they don't belong to the scene
			for (int i = 0; i < 4; i++) {
				std::getline(geomFile, lineStr);
			}
			// parse file based on given dimensions, will error if file does not match these
			// fills grid so that [0][0] is at bottom left corner of simulation
			for (int j = y - (borderCount + 1); j >= borderCount; j--) {
				std::getline(geomFile, lineStr);
				for (int i = borderCount; i < x - borderCount; i++) {
					switch (lineStr[i-borderCount]) {
					case 'f':
						grid.set(i, j, fixedZ, SimUtil::FLUID);
						break;
					case 's':
						grid.set(i, j, fixedZ, SimUtil::SOLID);
						break;
					case 'a':
						grid.set(i, j, fixedZ, SimUtil::AIR);
						break;
					}
				}
			}
			geomFile.close();
		}
	}

	void getDimensions(std::string FileName, int &width, int &height, int &depth, int &borderCount, int &emptyFrontCount, int &emptyBackCount, int &zDepth) {
		std::ifstream file(FileName);
		if (file.is_open()) {
			std::string lineStr;
			std::getline(file, lineStr);
			//thickness of the border arround the scene
			borderCount = std::stoi(lineStr);
			std::getline(file, lineStr);
			//number of empty z-layers in the front of the scene
			emptyFrontCount = std::stoi(lineStr);
			std::getline(file, lineStr);
			//number of empty z-layers in the back of the scene
			emptyBackCount = std::stoi(lineStr);
			std::getline(file, lineStr);
			//number of repetitions of the scene in z-layers
			zDepth = std::stoi(lineStr);
			depth = 2 * borderCount + emptyFrontCount + emptyBackCount + zDepth;
			std::getline(file, lineStr);
			int yCount = 1;
			width = 2 * borderCount + lineStr.length();
			while (std::getline(file, lineStr)) {
				yCount++;
			}
			//if file ends with a new line
			if (lineStr.length() < 2) {
				yCount--;
			}
			height = 2 * borderCount + yCount;
			file.close();
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