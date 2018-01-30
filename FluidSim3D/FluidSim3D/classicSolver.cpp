#include "classicSolver.h"
#include "SimUtil.h"
#include <iostream>

using namespace SimUtil;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

classicSolver::classicSolver(int gridWidth, int gridHeight, float dx, float dt, Mat2Di *label, Mat2Df *p, Mat2Df *u, Mat2Df *v)
{
	m_gridWidth = gridWidth;
	m_gridHeight = gridHeight;
	m_dx = dx;
	m_dt = dt;
	m_label = label;
	m_p = p;
	m_u = u;
	m_v = v;
}

//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------

classicSolver::~classicSolver()
{
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

/*
Solves for pressure using the current velocity field.
*/
void classicSolver::pressureSolve() {
	// initialize all grids to solve for pressure
	// using double for more accuracy
	Mat2Dd rhs{ m_gridWidth, m_gridHeight };
	constructRHS(rhs);
	Mat2Dd Adiag{ m_gridWidth, m_gridHeight };
	Mat2Dd Ax{ m_gridWidth, m_gridHeight };
	Mat2Dd Ay{ m_gridWidth, m_gridHeight };
	constructA(Adiag, Ax, Ay);
	Mat2Dd precon{ m_gridWidth, m_gridHeight };
	constructPrecon(precon, Adiag, Ax, Ay);

	// solve for pressure using PCG
	PCG(Adiag, Ax, Ay, rhs, precon);

	rhs.deleteGrid();
	Adiag.deleteGrid();
	Ax.deleteGrid();
	Ay.deleteGrid();
	precon.deleteGrid();
}

//----------------------------------------------------------------------
// Private Functions
//----------------------------------------------------------------------

/*
Sets up the right hand side of the system to solve for pressure. This is the negative
divergence at each cell center modified to account for the velocity of solids at boundaries.
Args:
rhs - the grid to use for the RHS
*/
void classicSolver::constructRHS(Mat2Dd &rhs) {
	// initialize to 0
	rhs.initValues(0.0);
	// calculate negative divergence
	double scale = 1.0f / m_dx;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label->get(i, j) == FLUID) {
				rhs.set(i,j, -scale * (m_u->get(i + 1, j) - m_u->get(i, j) + m_v->get(i, j + 1) - m_v->get(i, j)));
				// if it's on boundary must update to consider solid velocity
				// TODO create actual solid velocity grids, for right now just 0
				if (m_label->get(i - 1, j) == SOLID) {
					rhs.set(i, j, rhs.get(i,j) - scale * (m_u->get(i, j) - 0.0f)); //m_usolid[i][j]
				}
				if (m_label->get(i + 1, j) == SOLID) {
					rhs.set(i, j, rhs.get(i,j) + scale * (m_u->get(i + 1, j) - 0.0f)); //m_usolid[i+1][j]
				}
				if (m_label->get(i, j - 1) == SOLID) {
					rhs.set(i, j, rhs.get(i,j) - scale * (m_v->get(i, j) - 0.0f)); //m_vsolid[i][j]
				}
				if (m_label->get(i, j + 1) == SOLID) {
					rhs.set(i, j, rhs.get(i, j) + scale * (m_v->get(i, j + 1) - 0.0f)); //m_vsolid[i][j+1]
				}
			}
		}
	}
}

/*
Constructs the A matrix for the system to solve for pressure. This a sparse coefficient matrix
for the pressure terms, stored in 3 separate grids. If index i, j, k is not a fluid cell, then
it is 0.0 in all 3 grids that store the matix.
Args:
Adiag - grid to store the diagonal of the matrix in.
Ax - grid to store the coefficients for pressure in the (i+1) cell for each grid cell with x index i
Ay - grid to store the coefficients for pressure in the (j+1) cell for each grid cell with y index j
*/
void classicSolver::constructA(Mat2Dd &Adiag, Mat2Dd &Ax, Mat2Dd &Ay) {
	// clear to all zeros so can increment
	Adiag.initValues(0.0);
	Ax.initValues(0.0);
	Ay.initValues(0.0);

	// populate with coefficients for pressure unknowns
	double scale = m_dt / (FLUID_DENSITY * m_dx * m_dx);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label->get(i, j) == FLUID) {
				// handle negative x neighbor
				if (m_label->get(i - 1, j) == FLUID || m_label->get(i - 1, j) == AIR) {
					Adiag.set(i, j, Adiag.get(i,j) + scale);
				}
				// handle positive x neighbor
				if (m_label->get(i + 1, j) == FLUID) {
					Adiag.set(i, j, Adiag.get(i,j) + scale);
					Ax.set(i, j, -scale);
				}
				else if (m_label->get(i + 1, j) == AIR) {
					Adiag.set(i, j, Adiag.get(i, j) + scale);
				}
				// handle negative y neighbor
				if (m_label->get(i, j - 1) == FLUID || m_label->get(i, j - 1) == AIR) {
					Adiag.set(i, j, Adiag.get(i, j) + scale);
				}
				// handle positive y neighbor
				if (m_label->get(i, j + 1) == FLUID) {
					Adiag.set(i, j, Adiag.get(i, j) + scale);
					Ay.set(i, j, -scale);
				}
				else if (m_label->get(i, j + 1) == AIR) {
					Adiag.set(i, j, Adiag.get(i, j) + scale);
				}
			}
		}
	}
}

/*
Constructs the preconditioner used when performing the preconditioned conjugate gradient (PCG)
algorithm to solve for pressure.
Args:
precon - grid to store the preconditioner in
Adiag, Ax, Ay - the grids that make up the A coefficient matrix
*/
void classicSolver::constructPrecon(Mat2Dd &precon, Mat2Dd &Adiag, Mat2Dd &Ax, Mat2Dd &Ay) {
	precon.initValues(0.0);

	// tuning constant
	double tau = 0.97;
	// safety constant
	double sigma = 0.25;

	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label->get(i, j) == FLUID) {
				double Adiag_ij = Adiag.get(i, j);
				double Ax_im1j = 0.0;
				double Ax_ijm1 = 0.0;
				double Ay_ijm1 = 0.0;
				double Ay_im1j = 0.0;
				double precon_im1j = 0.0;
				double precon_ijm1 = 0.0;
				// want to stay at zero if off the grid
				// non-fluid entries in A are already 0
				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					if (m_label->get(i - 1, j) == FLUID) {
						Ax_im1j = Ax.get(i - 1, j);
						Ay_im1j = Ay.get(i - 1, j);
						precon_im1j = precon.get(i - 1, j);
					}
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					if (m_label->get(i, j - 1) == FLUID) {
						Ax_ijm1 = Ax.get(i, j - 1);
						Ay_ijm1 = Ay.get(i, j - 1);
						precon_ijm1 = precon.get(i, j - 1);
					}
				}

				double e = Adiag_ij - pow(Ax_im1j * precon_im1j, 2.0)
					- pow(Ay_ijm1 * precon_ijm1, 2.0)
					- tau * (
						Ax_im1j * Ay_im1j * pow(precon_im1j, 2.0)
						+ Ay_ijm1 * Ax_ijm1 * pow(precon_ijm1, 2.0)
						);

				if (e < (sigma * Adiag_ij)) {
					e = Adiag_ij;
				}

				precon.set(i, j, 1.0 / sqrt(e));
			}
		}
	}
}

/*
Performs the Modified Incomplete Cholesky Conjugate Gradient Level Zero
(preconditioned conjugate gradient, PCG) algorithm to solve the linear system
Ap = b for p. The result are placed on the pressure grid.
Args:
Adiag, Ax, Ay - the grids that make up the A coefficient matrix
b - the right hand side of the equation
precon - the preconditioner to use
*/
void classicSolver::PCG(Mat2Dd &Adiag, Mat2Dd &Ax, Mat2Dd &Ay, Mat2Dd &b, Mat2Dd &precon) {
	// reset pressure vec to 0
	m_p->initValues(0.0f);
	// residiual vector is b to start
	Mat2Dd r{ m_gridWidth, m_gridHeight };
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			r.set(i, j, b.get(i, j));
		}
	}

	// check if r = 0
	bool r0 = true;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (r.get(i, j) != 0) {
				r0 = false;
				break;
			}
		}
	}

	if (r0) {
		std::cout << "Did not run PCG b/c of 0 residual to start.\n";
		r.deleteGrid();
		return;
	}

	// auxiliary vector
	Mat2Dd z{ m_gridWidth, m_gridHeight };
	// search vector
	Mat2Dd s{ m_gridWidth, m_gridHeight };

	// initialize auxiliary and search
	applyPrecon(z, r, precon, Adiag, Ax, Ay);
	// search is set to auxiliary to start
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			s.set(i, j, z.get(i, j));
		}
	}

	double sigma = z.dot(r);

	// start the main iteration until tolerance is reach or max out iterations
	bool converged = false;
	for (int iters = 0; iters < PCG_MAX_ITERS; iters++) {
		applyA(z, s, Adiag, Ax, Ay);
		double alpha = sigma / z.dot(s);
		// update pressure and residual
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				m_p->set(i, j, m_p->get(i, j) + (alpha * s.get(i, j)));
				r.set(i, j, r.get(i, j) - (alpha * z.get(i, j)));
			}
		}
		// check if we're under the tolerance
		if (r.max() <= PCG_TOL) {
			//std::cout << "PCG converged after " << iters << " iterations.\n";
			converged = true;
			break;
		}
		// otherwise new auxiliary vector
		applyPrecon(z, r, precon, Adiag, Ax, Ay);
		double newSigma = z.dot(r);
		double beta = newSigma / sigma;
		// update search vector
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				s.set(i, j, z.get(i, j) + (beta * s.get(i, j)));
			}
		}
		// update sigma
		sigma = newSigma;
	}

	if (!converged) {
		std::cout << "PCG did not converge, stopped after " << PCG_MAX_ITERS << " iterations!\n";
	}

	r.deleteGrid();
	z.deleteGrid();
	s.deleteGrid();
}

/*
Applies the given preconditioner to the given vector and places the result in given "vector" z.
Args:
z - the vector (2D grid) to place the result in
r - the vector (2D grid) to multiply the preconditioner by
precon - the preconditioner
Adiag, Ax, Ay - the grids that make up the coefficient matrix
*/
void classicSolver::applyPrecon(Mat2Dd &z, Mat2Dd &r, Mat2Dd &precon, Mat2Dd &Adiag, Mat2Dd &Ax, Mat2Dd &Ay) {
	// first solve Lq = r
	Mat2Dd q{ m_gridWidth, m_gridHeight };
	q.initValues(0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label->get(i, j) == FLUID) {
				double Ax_im1j = 0.0;
				double Ay_ijm1 = 0.0;
				double precon_im1j = 0.0;
				double precon_ijm1 = 0.0;
				double q_im1j = 0.0;
				double q_ijm1 = 0.0;

				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					if (m_label->get(i - 1, j) == FLUID) {
						Ax_im1j = Ax.get(i - 1, j);
						precon_im1j = precon.get(i - 1, j);
						q_im1j = q.get(i - 1, j);
					}
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					if (m_label->get(i, j - 1) == FLUID) {
						Ay_ijm1 = Ay.get(i, j - 1);
						precon_ijm1 = precon.get(i, j - 1);
						q_ijm1 = q.get(i, j - 1);
					}
				}

				double t = r.get(i, j) - (Ax_im1j * precon_im1j * q_im1j)
					- (Ay_ijm1 * precon_ijm1 * q_ijm1);

				q.set(i, j, t * precon.get(i, j));
			}
		}
	}

	// now solve L^T z = q
	z.initValues(0.0);
	for (int i = m_gridWidth - 1; i >= 0; i--) {
		for (int j = m_gridHeight; j >= 0; j--) {
			if (m_label->get(i, j) == FLUID) {
				double Ax_ij = Ax.get(i, j);
				double Ay_ij = Ay.get(i, j);
				double precon_ij = precon.get(i, j);
				double z_ip1j = 0.0;
				double z_ijp1 = 0.0;

				if (i + 1 >= 0 && i + 1 < m_gridWidth) {
					if (m_label->get(i + 1, j) == FLUID) {
						z_ip1j = z.get(i + 1, j);
					}
				}
				if (j + 1 >= 0 && j + 1 < m_gridHeight) {
					if (m_label->get(i, j + 1) == FLUID) {
						z_ijp1 = z.get(i, j + 1);
					}
				}

				double t = q.get(i, j) - (Ax_ij * precon_ij * z_ip1j)
					- (Ay_ij * precon_ij * z_ijp1);

				z.set(i, j, t * precon_ij);
			}
		}
	}

	q.deleteGrid();
}

/*
Multiplies A by the given vector (a 2D grid).
Args:
z - the vector (2D grid) to place the result in
s - the vector (2D grid) to multiply A by
Adiag, Ax, Ay - the grids that make up the coefficient matrix A
*/
void classicSolver::applyA(Mat2Dd &z, Mat2Dd &s, Mat2Dd &Adiag, Mat2Dd &Ax, Mat2Dd &Ay) {
	z.initValues(0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label->get(i, j) == FLUID) {
				z.set(i, j, Adiag.get(i, j) * s.get(i, j)
					+ Ax.get(i, j) * s.get(i + 1, j)
					+ Ay.get(i, j) * s.get(i, j + 1));
				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					z.set(i, j, z.get(i, j) + Ax.get(i - 1, j) * s.get(i - 1, j));
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					z.set(i, j, z.get(i, j) + Ay.get(i, j - 1) * s.get(i, j - 1));
				}
			}
		}
	}
}

//----------------------------------------------------------------------
// Private Helper Functions
//----------------------------------------------------------------------

