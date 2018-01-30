#pragma once

#include "SimUtil.h"

class classicSolver
{
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
	SimUtil::Mat2Di *m_label;

	// pressure and velocity are held in a MAC grid so that
	// p(i, j, k) = p_i_j_k
	// u(i, j, k) = u_i-1/2_j_k
	// v(i, j, k) = v_i_j-1/2_k

	// grid of pressures, size (nx, ny)
	SimUtil::Mat2Df *m_p;
	// grid of vel x component, size (nx+1, ny)
	SimUtil::Mat2Df *m_u;
	// grid of vel y component, size (nx, ny+1)
	SimUtil::Mat2Df *m_v;


	// density of the fluid (kg/m^3)
	const float FLUID_DENSITY = 1000.0f;
	// error tolerance for PCG
	const float PCG_TOL = 0.000001f;
	// max iterations for PCG
	const int PCG_MAX_ITERS = 200;

	// simulation time step
	float m_dt;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------
	void constructRHS(SimUtil::Mat2Dd&);
	void constructA(SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&);
	void constructPrecon(SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&);
	void PCG(SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&);
	void applyPrecon(SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&);
	void applyA(SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&, SimUtil::Mat2Dd&);
public:
	classicSolver(int, int, float, float, SimUtil::Mat2Di*, SimUtil::Mat2Df*, SimUtil::Mat2Df*, SimUtil::Mat2Df*);
	~classicSolver();

	void pressureSolve();
};

