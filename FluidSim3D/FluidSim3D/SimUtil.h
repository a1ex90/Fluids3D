#ifndef SIM_UTIL_H
#define SIM_UTIL_H

#include <string>
#include <vector>
#include <glm/glm.hpp>

namespace SimUtil {
	class Mat2Di {
	public:
		Mat2Di();
		Mat2Di(int width, int height);
		void initValues(int val);
		void set(int width, int height, int val);
		int get(int width, int height);
		double dot(Mat2Di mat);
		int width();
		int height();
		int max();
		void deleteGrid();
		~Mat2Di();
	private:
		int m_width;
		int m_height;
		int* m_grid;
	};
	class Mat2Df {
	public:
		Mat2Df();
		Mat2Df(int width, int height);
		void initValues(float val);
		void set(int width, int height, float val);
		float get(int width, int height);
		double dot(Mat2Df mat);
		int width();
		int height();
		float max();
		void deleteGrid();
		~Mat2Df();
	private:
		int m_width;
		int m_height;
		float* m_grid;
	};
	class Mat2Dd {
	public:
		Mat2Dd();
		Mat2Dd(int width, int height);
		void initValues(double val);
		void set(int width, int height, double val);
		double get(int width, int height);
		double dot(Mat2Dd mat);
		int width();
		int height();
		double max();
		double *data();
		void deleteGrid();
		~Mat2Dd();
	private:
		int m_width;
		int m_height;
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

	struct Vec2 {
		float x, y;
		Vec2(float x, float y) : x(x), y(y) {}
	};

	struct Particle2D {
		Vec2 pos;
		Vec2 vel;
		Particle2D(Vec2 pos, Vec2 vel) : pos(pos), vel(vel) {}
	};

	struct MarchingTrianglesData {
		std::vector<glm::vec2> vertices;
		std::vector<int> indices;
		std::vector<float> opacities;
		MarchingTrianglesData(std::vector<glm::vec2> vert, std::vector<int> ind, std::vector<float> opa) : vertices(vert), indices(ind), opacities(opa) {}
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
	void readInGeom2D(int, int, std::string, SimUtil::Mat2Di&);

	/*
	Finds the physical location of the cell with index [x][y]
	in a grid where [0][0] is the center of the bottom left cell, based on the given dx. 
	Integer indices are treated in the center of cells while fractional indices may lie anywhere
	in a grid cell.
	Args:
	i - x index of cell
	j - y index of cell
	dx - single cell dimension
	Returns:
	Vec2 (x, y) containing physical location from bottom left corner of grid.
	*/
	Vec2 getGridCellPosition(float, float, float);
	/*
	Returns array of size 2 that contains the integer grid cell with index [i][j] at its center that contains the given position.
	Args:
	pos - a Vec2 (x,y) coordinate containing the position to use based on the origin at the bottom left of the simulation.
	dx - single cell dimension
	*/
	int* getGridCellIndex(Vec2 pos, float);

	// vec operations

	/*
	Calculates the sum of vec1 and vec2 (vec1 + vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec2 add(Vec2, Vec2);
	/*
	Calculates the difference of vec1 and vec2 (vec1 - vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec2 sub(Vec2, Vec2);
	/*
	Scales the given vector by a scalar and returns the new scaled vector. 
	Args
	vec1 - the vector
	scalar - the value to scale by
	*/
	Vec2 scale(Vec2, float);
	/*
	Calculates Euclidean norm of the given vector.
	Args
	vec - the vector to calculate norm of.
	*/
	float norm(Vec2);
}

#endif //SIM_UTIL_H

