#pragma once
#include <iostream>
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
	// nz
	int m_gridDepth;
	// distance between each grid cell
	float m_dx;
	// grid of cell labels, size (nx, ny)
	SimUtil::Mat3Di *m_label;

	// pressure and velocity are held in a MAC grid so that
	// p(i, j, k) = p_i_j_k
	// u(i, j, k) = u_i-1/2_j_k
	// v(i, j, k) = v_i_j-1/2_k
	// w(i, j, k) = w_i_j_k-1/2

	// grid of pressures, size (nx, ny)
	SimUtil::Mat3Df *m_p;
	// grid of vel x component, size (nx+1, ny, nz)
	SimUtil::Mat3Df *m_u;
	// grid of vel y component, size (nx, ny+1, nz)
	SimUtil::Mat3Df *m_v;
	// grid of vel z component, size (nx, ny, nz+1)
	SimUtil::Mat3Df *m_w;


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
	void constructRHS(SimUtil::Mat3Dd&);
	void constructA(SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&);
	void constructPrecon(SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&);
	void PCG(SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&);
	void applyPrecon(SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&, SimUtil::Mat3Dd&);
	void applyA(SimUtil::Mat3Dd &z, SimUtil::Mat3Dd &s, SimUtil::Mat3Dd &Adiag, SimUtil::Mat3Dd &Ax, SimUtil::Mat3Dd &Ay, SimUtil::Mat3Dd &Az);
public:
	classicSolver(int, int, int, float, float, SimUtil::Mat3Di*, SimUtil::Mat3Df*, SimUtil::Mat3Df*, SimUtil::Mat3Df*, SimUtil::Mat3Df*);
	~classicSolver();

	void pressureSolve();
};

