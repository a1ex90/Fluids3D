#ifndef SIM_UTIL_H
#define SIM_UTIL_H

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

namespace SimUtil {
	class Mat3Di {
	public:
		Mat3Di();
		Mat3Di(int width, int height, int depth);
		void initValues(int val);
		void set(int width, int height, int depth, int val);
		int get(int width, int height, int depth);
		double dot(Mat3Di mat);
		int width();
		int height();
		int depth();
		int max();
		void deleteGrid();
		~Mat3Di();
	private:
		int m_width;
		int m_height;
		int m_depth;
		int* m_grid;
	};
	class Mat3Df {
	public:
		Mat3Df();
		Mat3Df(int width, int height, int depth);
		void initValues(float val);
		void set(int width, int height, int depth, float val);
		float get(int width, int height, int depth);
		double dot(Mat3Df mat);
		int width();
		int height();
		int depth();
		float max();
		void deleteGrid();
		~Mat3Df();
	private:
		int m_width;
		int m_height;
		int m_depth;
		float* m_grid;
	};
	class Mat3Dd {
	public:
		Mat3Dd();
		Mat3Dd(int width, int height, int depth);
		void initValues(double val);
		void set(int width, int height, int depth, double val);
		double get(int width, int height, int depth);
		double dot(Mat3Dd mat);
		int width();
		int height();
		int depth();
		double max();
		double *data();
		void deleteGrid();
		~Mat3Dd();
	private:
		int m_width;
		int m_height;
		int m_depth;
		double* m_grid;
	};

	//----------------------------------------------------------------------
	// Constants
	//----------------------------------------------------------------------

	const int SOLID = 0;
	const int FLUID = 1;
	const int AIR = 2;

	//----------------------------------------------------------------------
	// Data Structures
	//----------------------------------------------------------------------

	struct Vec3 {
		float x, y, z;
		Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
	};

	struct Particle3D {
		Vec3 pos;
		Vec3 vel;
		Particle3D(Vec3 pos, Vec3 vel) : pos(pos), vel(vel) {}
	};

	struct MarchingTrianglesData {
		std::vector<glm::vec2> vertices;
		std::vector<int> indices;
		std::vector<float> opacities;
		MarchingTrianglesData(std::vector<glm::vec2> vert, std::vector<int> ind, std::vector<float> opa) : vertices(vert), indices(ind), opacities(opa) {}
	};

	struct Mesh3D {
		std::vector<glm::vec3> vertices;
		std::vector<glm::vec3> normals;
		std::vector<int> indices;
		Mesh3D(std::vector<glm::vec3> vert, std::vector<glm::vec3> norm, std::vector<int> ind) : vertices(vert), normals(norm), indices(ind) {}
	};

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	/*
	Builds initial grid of dimensions (x, y) that contains the initial
	geometry for the system to simulate. Cell [0][0] in the grid is at the bottom
	left corner of the input geometry, so it's treated as if the input grid was initialized
	using initGrid2D. x is positive right, y is positive up. 
	It reads the initial geometry from the specified input file parameter.
	Args:
	x, y - grid dimensions in number of cells
	geomFile - the file containing the geometry
	grid - the 2D array to put the initial grid in
	*/
	void readInGeom2D(int x, int y, int fixedZ, int borderCount, std::string geomFileName, Mat3Di &grid);
	void readInGeom3D(int x, int y, int z, int borderCount, int emptyFrontCount, int emptyBackCount, int zDepth, std::string geomFileName, Mat3Di &grid);

	void getDimensions(std::string FileName, int &width, int &height, int &depth, int &borderCount, int &emptyFrontCount, int &emptyBackCount, int &zDepth);

	/*
	Finds the physical location of the cell with index [x][y][z]
	in a grid where [0][0][0] is the center of the bottom left cell, based on the given dx. 
	Integer indices are treated in the center of cells while fractional indices may lie anywhere
	in a grid cell.
	Args:
	i - x index of cell
	j - y index of cell
	k - z index of cell
	dx - single cell dimension
	Returns:
	Vec3 (x, y, z) containing physical location from bottom left corner of grid.
	*/
	Vec3 getGridCellPosition(float, float, float, float);
	/*
	Returns array of size 2 that contains the integer grid cell with index [i][j] at its center that contains the given position.
	Args:
	pos - a Vec2 (x,y) coordinate containing the position to use based on the origin at the bottom left of the simulation.
	dx - single cell dimension
	*/
	int* getGridCellIndex(Vec3 pos, float);

	// vec operations

	/*
	Calculates the sum of vec1 and vec2 (vec1 + vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec3 add(Vec3, Vec3);
	/*
	Calculates the difference of vec1 and vec2 (vec1 - vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec3 sub(Vec3, Vec3);
	/*
	Scales the given vector by a scalar and returns the new scaled vector. 
	Args
	vec1 - the vector
	scalar - the value to scale by
	*/
	Vec3 scale(Vec3, float);
	/*
	Calculates Euclidean norm of the given vector.
	Args
	vec - the vector to calculate norm of.
	*/
	float norm(Vec3);
}

#endif //SIM_UTIL_H

