#include "classicSolver.h"

using namespace SimUtil;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

classicSolver::classicSolver(int gridWidth, int gridHeight, int gridDepth, float dx, float dt, Mat3Di *label, Mat3Df *p, Mat3Df *u, Mat3Df *v, Mat3Df *w)
{
	m_gridWidth = gridWidth;
	m_gridHeight = gridHeight;
	m_gridDepth = gridDepth;
	m_dx = dx;
	m_dt = dt;
	m_label = label;
	m_p = p;
	m_u = u;
	m_v = v;
	m_w = w;
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
	Mat3Dd rhs{ m_gridWidth, m_gridHeight, m_gridDepth };
	constructRHS(rhs);
	Mat3Dd Adiag{ m_gridWidth, m_gridHeight, m_gridDepth };
	Mat3Dd Ax{ m_gridWidth, m_gridHeight, m_gridDepth };
	Mat3Dd Ay{ m_gridWidth, m_gridHeight, m_gridDepth };
	Mat3Dd Az{ m_gridWidth, m_gridHeight, m_gridDepth };
	constructA(Adiag, Ax, Ay, Az);
	Mat3Dd precon{ m_gridWidth, m_gridHeight, m_gridDepth };
	constructPrecon(precon, Adiag, Ax, Ay, Az);

	// solve for pressure using PCG
	PCG(Adiag, Ax, Ay, Az, rhs, precon);

	rhs.deleteGrid();
	Adiag.deleteGrid();
	Ax.deleteGrid();
	Ay.deleteGrid();
	Az.deleteGrid();
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
void classicSolver::constructRHS(Mat3Dd &rhs) {
	// initialize to 0
	rhs.initValues(0.0);
	// calculate negative divergence
	double scale = 1.0f / m_dx;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label->get(i, j, k) == FLUID) {
					rhs.set(i, j, k, -scale * (m_u->get(i + 1, j, k) - m_u->get(i, j, k) + m_v->get(i, j + 1, k) - m_v->get(i, j, k) + m_w->get(i, j, k + 1) - m_w->get(i, j, k)));
					// if it's on boundary must update to consider solid velocity
					// TODO create actual solid velocity grids, for right now just 0
					if (m_label->get(i - 1, j, k) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) - scale * (m_u->get(i, j, k) - 0.0f)); //m_usolid[i][j]
					}
					if (m_label->get(i + 1, j, k) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) + scale * (m_u->get(i + 1, j, k) - 0.0f)); //m_usolid[i+1][j]
					}
					if (m_label->get(i, j - 1, k) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) - scale * (m_v->get(i, j, k) - 0.0f)); //m_vsolid[i][j]
					}
					if (m_label->get(i, j + 1, k) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) + scale * (m_v->get(i, j + 1, k) - 0.0f)); //m_vsolid[i][j+1]
					}

					if (m_label->get(i, j, k - 1) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) - scale * (m_w->get(i, j, k) - 0.0f)); //m_wsolid[i][j]
					}
					if (m_label->get(i, j, k + 1) == SOLID) {
						rhs.set(i, j, k, rhs.get(i, j, k) + scale * (m_w->get(i, j, k + 1) - 0.0f)); //m_wsolid[i][j+1]
					}
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
Az - grid to store the coefficients for pressure in the (k+1) cell for each grid cell with z index k
*/
void classicSolver::constructA(Mat3Dd &Adiag, Mat3Dd &Ax, Mat3Dd &Ay, Mat3Dd &Az) {
	// clear to all zeros so can increment
	Adiag.initValues(0.0);
	Ax.initValues(0.0);
	Ay.initValues(0.0);
	Az.initValues(0.0);

	// populate with coefficients for pressure unknowns
	double scale = m_dt / (FLUID_DENSITY * m_dx * m_dx);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label->get(i, j, k) == FLUID) {
					// handle negative x neighbor
					if (m_label->get(i - 1, j, k) == FLUID || m_label->get(i - 1, j, k) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
					// handle positive x neighbor
					if (m_label->get(i + 1, j, k) == FLUID) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
						Ax.set(i, j, k, -scale);
					}
					else if (m_label->get(i + 1, j, k) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
					// handle negative y neighbor
					if (m_label->get(i, j - 1, k) == FLUID || m_label->get(i, j - 1, k) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
					// handle positive y neighbor
					if (m_label->get(i, j + 1, k) == FLUID) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
						Ay.set(i, j, k, -scale);
					}
					else if (m_label->get(i, j + 1, k) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
					// handle negative z neighbor
					if (m_label->get(i, j, k - 1) == FLUID || m_label->get(i, j, k - 1) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
					// handle positive z neighbor
					if (m_label->get(i, j, k + 1) == FLUID) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
						Az.set(i, j, k, -scale);
					}
					else if (m_label->get(i, j, k + 1) == AIR) {
						Adiag.set(i, j, k, Adiag.get(i, j, k) + scale);
					}
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
Adiag, Ax, Ay, Az - the grids that make up the A coefficient matrix
*/
void classicSolver::constructPrecon(Mat3Dd &precon, Mat3Dd &Adiag, Mat3Dd &Ax, Mat3Dd &Ay, Mat3Dd &Az) {
	precon.initValues(0.0);

	// tuning constant
	double tau = 0.97;
	// safety constant
	double sigma = 0.25;

	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label->get(i, j, k) == FLUID) {
					double Adiag_ijk = Adiag.get(i, j, k);
					double Ax_im1jk = 0.0;
					double Ax_ijm1k = 0.0;
					double Ax_ijkm1 = 0.0;
					double Ay_ijm1k = 0.0;
					double Ay_im1jk = 0.0;
					double Ay_ijkm1 = 0.0;
					double Az_im1jk = 0.0;
					double Az_ijm1k = 0.0;
					double Az_ijkm1 = 0.0;
					double precon_im1jk = 0.0;
					double precon_ijm1k = 0.0;
					double precon_ijkm1 = 0.0;
					// want to stay at zero if off the grid
					// non-fluid entries in A are already 0
					if (i - 1 >= 0 && i - 1 < m_gridWidth) {
						if (m_label->get(i - 1, j, k) == FLUID) {
							Ax_im1jk = Ax.get(i - 1, j, k);
							Ay_im1jk = Ay.get(i - 1, j, k);
							Az_im1jk = Az.get(i - 1, j, k);
							precon_im1jk = precon.get(i - 1, j, k);
						}
					}
					if (j - 1 >= 0 && j - 1 < m_gridHeight) {
						if (m_label->get(i, j - 1, k) == FLUID) {
							Ax_ijm1k = Ax.get(i, j - 1, k);
							Ay_ijm1k = Ay.get(i, j - 1, k);
							Az_ijm1k = Az.get(i, j - 1, k);
							precon_ijm1k = precon.get(i, j - 1, k);
						}
					}
					if (k - 1 >= 0 && k - 1 < m_gridDepth) {
						if (m_label->get(i, j, k - 1) == FLUID) {
							Ax_ijkm1 = Ax.get(i, j, k - 1);
							Ay_ijkm1 = Ay.get(i, j, k - 1);
							Az_ijkm1 = Az.get(i, j, k - 1);
							precon_ijkm1 = precon.get(i, j, k - 1);
						}
					}

					double e = Adiag_ijk - (Ax_im1jk * precon_im1jk)*(Ax_im1jk * precon_im1jk)
						- (Ay_ijm1k * precon_ijm1k) * (Ay_ijm1k * precon_ijm1k)
						- (Az_ijkm1 * precon_ijkm1) * (Az_ijkm1 * precon_ijkm1)
						- tau * (
							Ax_im1jk * (Ay_im1jk + Az_im1jk) * precon_im1jk * precon_im1jk
							+ Ay_ijm1k * (Ax_ijm1k + Az_ijm1k) * precon_ijm1k * precon_ijm1k
							+ Az_ijkm1 * (Ax_ijkm1 + Ay_ijkm1) * precon_ijkm1 * precon_ijkm1
							);

					if (e < (sigma * Adiag_ijk)) {
						e = Adiag_ijk;
					}

					precon.set(i, j, k, 1.0 / sqrt(e));
				}
			}
		}
	}
}

/*
Performs the Modified Incomplete Cholesky Conjugate Gradient Level Zero
(preconditioned conjugate gradient, PCG) algorithm to solve the linear system
Ap = b for p. The result are placed on the pressure grid.
Args:
Adiag, Ax, Ay, Az - the grids that make up the A coefficient matrix
b - the right hand side of the equation
precon - the preconditioner to use
*/
void classicSolver::PCG(Mat3Dd &Adiag, Mat3Dd &Ax, Mat3Dd &Ay, Mat3Dd &Az, Mat3Dd &b, Mat3Dd &precon) {
	// reset pressure vec to 0
	m_p->initValues(0.0f);
	// residiual vector is b to start
	Mat3Dd r{ m_gridWidth, m_gridHeight, m_gridDepth };
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				r.set(i, j, k, b.get(i, j, k));
			}
		}
	}

	// check if r = 0
	bool r0 = true;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (r.get(i, j, k) != 0) {
					r0 = false;
					break;
				}
			}
		}
	}

	if (r0) {
		std::cout << "Did not run PCG b/c of 0 residual to start.\n";
		r.deleteGrid();
		return;
	}

	// auxiliary vector
	Mat3Dd z{ m_gridWidth, m_gridHeight, m_gridDepth };
	// search vector
	Mat3Dd s{ m_gridWidth, m_gridHeight, m_gridDepth };

	// initialize auxiliary and search
	applyPrecon(z, r, precon, Adiag, Ax, Ay, Az);
	// search is set to auxiliary to start
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				s.set(i, j, k, z.get(i, j, k));
			}
		}
	}

	double sigma = z.dot(r);

	// start the main iteration until tolerance is reach or max out iterations
	bool converged = false;
	for (int iters = 0; iters < PCG_MAX_ITERS; iters++) {
		applyA(z, s, Adiag, Ax, Ay, Az);
		double alpha = sigma / z.dot(s);
		// update pressure and residual
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				for (int k = 0; k < m_gridDepth; k++) {
					m_p->set(i, j, k, m_p->get(i, j, k) + (alpha * s.get(i, j, k)));
					r.set(i, j, k, r.get(i, j, k) - (alpha * z.get(i, j, k)));
				}
			}
		}
		// check if we're under the tolerance
		if (r.max() <= PCG_TOL) {
			//std::cout << "PCG converged after " << iters << " iterations.\n";
			converged = true;
			break;
		}
		// otherwise new auxiliary vector
		applyPrecon(z, r, precon, Adiag, Ax, Ay, Az);
		double newSigma = z.dot(r);
		double beta = newSigma / sigma;
		// update search vector
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				for (int k = 0; k < m_gridDepth; k++) {
					s.set(i, j, k, z.get(i, j, k) + (beta * s.get(i, j, k)));
				}
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
z - the vector (3D grid) to place the result in
r - the vector (3D grid) to multiply the preconditioner by
precon - the preconditioner
Adiag, Ax, Ay, Az - the grids that make up the coefficient matrix
*/
void classicSolver::applyPrecon(Mat3Dd &z, Mat3Dd &r, Mat3Dd &precon, Mat3Dd &Adiag, Mat3Dd &Ax, Mat3Dd &Ay, Mat3Dd &Az) {
	// first solve Lq = r
	Mat3Dd q{ m_gridWidth, m_gridHeight, m_gridDepth };
	q.initValues(0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label->get(i, j, k) == FLUID) {
					double Ax_im1jk = 0.0;
					double Ay_ijm1k = 0.0;
					double Az_ijkm1 = 0.0;
					double precon_im1jk = 0.0;
					double precon_ijm1k = 0.0;
					double precon_ijkm1 = 0.0;
					double q_im1jk = 0.0;
					double q_ijm1k = 0.0;
					double q_ijkm1 = 0.0;

					if (i - 1 >= 0 && i - 1 < m_gridWidth) {
						if (m_label->get(i - 1, j, k) == FLUID) {
							Ax_im1jk = Ax.get(i - 1, j, k);
							precon_im1jk = precon.get(i - 1, j, k);
							q_im1jk = q.get(i - 1, j, k);
						}
					}
					if (j - 1 >= 0 && j - 1 < m_gridHeight) {
						if (m_label->get(i, j - 1, k) == FLUID) {
							Ay_ijm1k = Ay.get(i, j - 1, k);
							precon_ijm1k = precon.get(i, j - 1, k);
							q_ijm1k = q.get(i, j - 1, k);
						}
					}
					if (k - 1 >= 0 && k - 1 < m_gridDepth) {
						if (m_label->get(i, j, k - 1) == FLUID) {
							Az_ijkm1 = Az.get(i, j, k - 1);
							precon_ijkm1 = precon.get(i, j, k - 1);
							q_ijkm1 = q.get(i, j, k - 1);
						}
					}

					double t = r.get(i, j, k) 
						- (Ax_im1jk * precon_im1jk * q_im1jk)
						- (Ay_ijm1k * precon_ijm1k * q_ijm1k)
						- (Az_ijkm1 * precon_ijkm1 * q_ijkm1);

					q.set(i, j, k, t * precon.get(i, j, k));
				}
			}
		}
	}

	// now solve L^T z = q
	z.initValues(0.0);
	for (int i = m_gridWidth - 1; i >= 0; i--) {
		for (int j = m_gridHeight - 1; j >= 0; j--) {
			for (int k = m_gridDepth - 1; k >= 0; k--) {
				if (m_label->get(i, j, k) == FLUID) {
					double Ax_ijk = Ax.get(i, j, k);
					double Ay_ijk = Ay.get(i, j, k);
					double Az_ijk = Az.get(i, j, k);
					double precon_ijk = precon.get(i, j, k);
					double z_ip1jk = 0.0;
					double z_ijp1k = 0.0;
					double z_ijkp1 = 0.0;

					if (i + 1 >= 0 && i + 1 < m_gridWidth) {
						if (m_label->get(i + 1, j, k) == FLUID) {
							z_ip1jk = z.get(i + 1, j, k);
						}
					}
					if (j + 1 >= 0 && j + 1 < m_gridHeight) {
						if (m_label->get(i, j + 1, k) == FLUID) {
							z_ijp1k = z.get(i, j + 1, k);
						}
					}
					if (k + 1 >= 0 && k + 1 < m_gridDepth) {
						if (m_label->get(i, j, k + 1) == FLUID) {
							z_ijkp1 = z.get(i, j, k + 1);
						}
					}

					double t = q.get(i, j, k) 
						- (Ax_ijk * precon_ijk * z_ip1jk)
						- (Ay_ijk * precon_ijk * z_ijp1k)
						- (Az_ijk * precon_ijk * z_ijkp1);

					z.set(i, j, k, t * precon_ijk);
				}
			}
		}
	}

	q.deleteGrid();
}

/*
Multiplies A by the given vector (a 2D grid).
Args:
z - the vector (3D grid) to place the result in
s - the vector (3D grid) to multiply A by
Adiag, Ax, Ay, Az - the grids that make up the coefficient matrix A
*/
void classicSolver::applyA(Mat3Dd &z, Mat3Dd &s, Mat3Dd &Adiag, Mat3Dd &Ax, Mat3Dd &Ay, Mat3Dd &Az) {
	z.initValues(0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			for (int k = 0; k < m_gridDepth; k++) {
				if (m_label->get(i, j, k) == FLUID) {
					z.set(i, j, k, Adiag.get(i, j, k) * s.get(i, j, k)
						+ Ax.get(i, j, k) * s.get(i + 1, j, k)
						+ Ay.get(i, j, k) * s.get(i, j + 1, k)
						+ Az.get(i, j, k) * s.get(i, j, k + 1));
					if (i - 1 >= 0 && i - 1 < m_gridWidth) {
						z.set(i, j, k, z.get(i, j, k) + Ax.get(i - 1, j, k) * s.get(i - 1, j, k));
					}
					if (j - 1 >= 0 && j - 1 < m_gridHeight) {
						z.set(i, j, k, z.get(i, j, k) + Ay.get(i, j - 1, k) * s.get(i, j - 1, k));
					}
					if (k - 1 >= 0 && k - 1 < m_gridDepth) {
						z.set(i, j, k, z.get(i, j, k) + Az.get(i, j, k - 1) * s.get(i, j, k - 1));
					}
				}
			}
		}
	}
}

//----------------------------------------------------------------------
// Private Helper Functions
//----------------------------------------------------------------------

