#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <iterator>
#include <fstream>

#include <data_type.h>

/*! The spline class that generates and stores a minimal spline trajectory. 

	MinimisingPoly uses Algorithm 1 to solve (1) to generate the spline trajectory.

	\tparam s The largest degree of continuity enforced by the spline.
*/
template <int s>
class MinimisingPoly {
    int K;
    double *T;
	// double *B;
	std::map<int, double> *B;
    VectorP<s> *P;
	VectorP<s> *d;
	MatrixBlock<s> *Hi;

	bool save_mem;
    bool debug;
public: 
	/*! Default constructor.
	
	\param number_segments Number of segments in the spline.
	\param save_mem_Data A flag (set false) to save intermediate matrix calculations.
	\param debug_data A flag (set true) for some debug output.
	*/
    MinimisingPoly(int number_segments, bool save_mem_data, bool debug_data) {
		if (debug_data) {
			std::cout << "Initialising object.\n";
		}
		K = number_segments;
		save_mem = save_mem_data;
		debug = debug_data;

		P = new VectorP<s> [K];
		d = new VectorP<s> [K];
		if (!save_mem) {
			Hi = new MatrixBlock<s> [K];
		}
	}

	MinimisingPoly(int number_segments) {
		debug = false;
		K = number_segments;

		P = new VectorP<s> [K];
		d = new VectorP<s> [K];

		save_mem = false;
		Hi = new MatrixBlock<s> [K];
	}

	~MinimisingPoly() {
		delete[] P;
		delete[] d;
		if (!save_mem) {
			delete[] Hi;
		}
	}

	/*! Solves (1) using Algorithm 1, storing coefficients of k segments of order 2*s at *P. 
	*/
    void solve_minimisation() {
		try {
			if (debug) {
				std::cout << "Allocating memory for matrices used in calculations.\n";
			}
			if (save_mem) {
				Hi = new MatrixBlock<s> [K];
			}
			MatrixSubBCD<s> *Ri = new MatrixSubBCD<s> [K];

			gen_block_matrices(Hi, Ri);

			partitioning(Hi, Ri);

			forward_updates(Hi, Ri);

			backwards_solution(Ri);

			solve_for_coefficients(Hi, Ri);

			if (save_mem) {
				delete[] Hi;
			}
			delete[] Ri;
		}
		catch (const std::bad_alloc& ba) {
			std::cerr << "bad_alloc caught: " << ba.what() << '\n';	
		}
	}

	/*! Creates the matrices and vectors that are the LHS of (3a) and (3b), respectively. 
	
		\param Hi Array of structs containing the output matrices and vectors.
		\param K Number of segments in spline.
		\param T Array of segment times.
		\param B Array of fixed derivative values.
	*/
	void gen_block_matrices(MatrixBlock<s> *Hi, MatrixSubBCD<s> *Ri) {
		MatrixQ<s> Q_temp;
		Q_temp = gen_Q_matrix(1, -1);
		for (int i = 0; i < K; i++) {
			// Q_temp = gen_Q_matrix(T[i + 1], T[i]);
			Hi[i].Q = Q_temp;
			MatrixA<s> A_temp = gen_Gamma_A_matrix(T[i + 1], T[i]);
			Hi[i].A = A_temp;
			MatrixA<s> A_temp_inv = A_temp.partialPivLu().inverse();
			Hi[i].H = (A_temp_inv.transpose()) * Q_temp * A_temp_inv;

			Ri[i].b = VectorB<s>::Zero();
			std::map<int, double>::iterator itr;
			for (itr = B[i].begin(); itr != B[i].end(); itr++) { 
				Ri[i].b(itr->first) = itr->second;
			}

			Hi[i].g = Hi[i].H * Ri[i].b;
			Ri[i].Z = BigMatrixPartition<s>::Identity(2*s,2*s);
		}
	}

	void partitioning(MatrixBlock<s> *Hi, MatrixSubBCD<s> *Ri) {
		for (int i = 0; i < K; i++) {
			fill_block(Ri[i], Hi[i].H, Hi[i].g);
			std::map<int, double>::reverse_iterator ritr;
			for (ritr = B[i].rbegin(); ritr != B[i].rend(); ritr++) { 
				reduce_block(Ri[i], ritr->first);
			}
		}
	}

	void partition_block(MatrixSubBCD<s> &R, MatrixA<s> &H, VectorB<s> &g) {
		fill_block(R, H, g);
		reduce_block(R, s);
		reduce_block(R, 0);
	}

	/*! Recursively updates the submatrices of (3a) and (3b) as the first for loop in Algorithm 1. 
	
		\param Ri Array of structs containing the partioned matrices and vectors.
		\param Hi Array of structs containing the unpartioned matrices and vectors.
		\param K Number of segments in spline.
	*/
	void forward_updates(MatrixBlock<s> *Hi, MatrixSubBCD<s> *Ri) {
		MatrixPartition<s> B_old;
		MatrixPartition<s> C_old;
		MatrixPartition<s> D_old;
		VectorPartition<s> g0_old;
		VectorPartition<s> gT_old;
		// MatrixPartition<s> B_old = MatrixPartition<s>::Zero(s-1, s-1);
		// MatrixPartition<s> C_old = MatrixPartition<s>::Zero(s-1, s-1);
		// MatrixPartition<s> D_old = MatrixPartition<s>::Zero(s-1, s-1);
		// VectorPartition<s> g0_old = VectorPartition<s>::Zero(s-1, s-1);
		// VectorPartition<s> gT_old = VectorPartition<s>::Zero(s-1, s-1);

		for (int i = 0; i < K; i++) {
			if (i == 0) {
				B_old = Ri[0].B; C_old = Ri[0].C; D_old = Ri[0].D;
				g0_old = Ri[0].g0; gT_old = Ri[0].gT;
			} else {
				MatrixPartition<s> B_inv_temp = B_old.partialPivLu().inverse();
				Ri[i].B = Ri[i].B + D_old - (C_old.transpose() * B_inv_temp * C_old);
				Ri[i].g0 = Ri[i].g0 + gT_old - (C_old.transpose() * B_inv_temp * g0_old);

				B_old = Ri[i].B; C_old = Ri[i].C; D_old = Ri[i].D;
				g0_old = Ri[i].g0; gT_old = Ri[i].gT;
			}
		}
	}

	/*! Recursively solves for the derivatives as the second for loop in Algorithm 1.
	
		\param Ri Array of structs containing the partioned matrices and vectors.
		\param f0 Array of vectors of the derivatives at the start of each segment.
		\param fT Array of vectors of the derivatives at the end of each segment.
		\param K Number of segments in spline.
	*/
	void backwards_solution(MatrixSubBCD<s> *Ri) {
		for (int i = K - 1; i >= 0; i--) {
			if (i == K - 1) {
				MatrixPartition<s> B_inv_temp = Ri[K - 1].B.partialPivLu().inverse();
				MatrixPartition<s> temp = Ri[K - 1].D - (Ri[K - 1].C.transpose() * B_inv_temp * Ri[K - 1].C);
				VectorPartition<s> v_temp = Ri[K - 1].gT - (Ri[K - 1].C.transpose() * B_inv_temp * Ri[K - 1].g0);
				Ri[K - 1].fT = temp.partialPivLu().solve(v_temp);
			} else {
				Ri[i].fT = Ri[i + 1].f0;
			}
			Ri[i].f0 = Ri[i].B.partialPivLu().solve(Ri[i].g0 - (Ri[i].C * Ri[i].fT));
		}
	}
	
	/*! Solves for polynomial coefficients by solving (1b).
	
		\param P Array of vectors of polynomial coefficients for each segment. 
		\param d Array of vectors of derivatives at beginning and end of each segment.
		\param H Array of structs containing matrices of Program (1).
		\param f0 Array of vectors of the derivatives at the start of each segment.
		\param fT Array of vectors of the derivatives at the end of each segment.
		\param K Number of segments in spline.
	*/
	void solve_for_coefficients(MatrixBlock<s> *H, MatrixSubBCD<s> *Ri) {
		for (int i = 0; i < K; i++) {
			BigMatrixPartition<s> f(Ri[i].f0.rows() + Ri[i].fT.rows(), 1);
			f.block(0, 0, Ri[i].f0.rows(), 1) = Ri[i].f0;
			f.block(Ri[i].f0.rows(), 0, Ri[i].fT.rows(), 1) = Ri[i].fT;
			VectorB<s> d_i = Ri[i].b + (Ri[i].Z * f);
			d[i] = d_i;
			P[i] = H[i].A.partialPivLu().solve(d_i);
		}
	}

	/*! Creates matrix A in Program (1).
	
		\param t1 Time at end of corresponding segment.
		\param t0 Time at beginning of corresponding segment.
	*/
	MatrixA<s> gen_Gamma_A_matrix(double t1, double t0) {
		MatrixA<s> A = MatrixA<s>::Zero();
		for (int h = 0; h < 2*s; h = h + s) {
			for (int i = 0; i < s; i++) {
				for (int j = i; j < 2*s; j++) {
					double coeff = 1;
					for (int k = i; k > 0; k--) {
						coeff = coeff * (j - k + 1);
					}
					if (h == 0) {
						A(i + h, j) = coeff * pow(-1, j - i) * pow(2.0 / (t1 - t0), i);
						// A(i + h, j) = coeff * pow(t0, j - i);
					}
					else if (h == s) {
						A(i + h, j) = coeff * pow(1, j - i) * pow(2.0 / (t1 - t0), i);
						// A(i + h, j) = coeff * pow(t1, j - i);
					}
				}
			}
		}
		return A;
	}

	/*! Creates matrix Q in Program (1).
	
		\param t1 Time at end of corresponding segment.
		\param t0 Time at beginning of corresponding segment.
	*/
	MatrixQ<s>gen_Q_matrix(double t1, double t0) {
		MatrixQ<s> Q = MatrixQ<s>::Zero();
		MatrixQ<s> Q_test1= MatrixQ<s>::Zero();
		VectorP<s> a_test1 = VectorP<s>::Zero();
		double t = t1;
		// double t = 1.0;
		// take outer product
		for (int i = int(s - 1); i < 2*s; i++) {
			double coeff = 1;
			for (int j = i; j > i - (s - 1); j--) {
				coeff = coeff * j;
			}
			a_test1(i) = coeff * pow(t, i - 1.0);
		}
		for (int i = int(s - 1); i < 2*s; i++) {
			for (int j = int(s - 1); j < 2*s; j++) {
				Q_test1(i, j) = a_test1(i) * a_test1(j) * (t / ((i - s + 1) + (j - s + 1) + 1.0));
			}
		}
		MatrixQ<s> Q_test_1 = MatrixQ<s>::Zero();
		VectorP<s> a_test_1 = VectorP<s>::Zero();
		// t = -1.0;
		t = t0;
		// take outer product
		for (int i = int(s - 1); i < 2*s; i++) {
			double coeff = 1;
			for (int j = i; j > i - (s - 1); j--) {
				coeff = coeff * j;
			}
			a_test_1(i) = coeff * pow(t, i - 1.0);
		}
		for (int i = int(s - 1); i < 2*s; i++) {
			for (int j = int(s - 1); j < 2*s; j++) {
				Q_test_1(i, j) = a_test_1(i) * a_test_1(j) * (t / ((i - s + 1) + (j - s + 1) + 1.0));
			}
		}
		Q = Q_test1 - Q_test_1;
		return Q;
	}

	/*! Partions the submatrices as (3a) and (3b).
	
		\param R Struct containing the submatrices.
		\param H LHS of (3a).
		\param g LHS of (3b).
	*/

	void fill_block(MatrixSubBCD<s> &R, MatrixA<s> &H, VectorB<s> &g) {
		R.B = H.template block<s,s>(0, 0);
		R.C = H.template block<s,s>(0, s);
		R.D = H.template block<s,s>(s, s);
		R.g0 = (-1) * g.template block<s,1>(0, 0);
		R.gT = (-1) * g.template block<s,1>(s, 0);
	}

	void reduce_block(MatrixSubBCD<s> &R, int del) {
		if ((0 <= del) && (del < s)) {
			remove_col(R.B, del);
			remove_row(R.B, del);
			remove_row(R.C, del);
			remove_row(R.g0, del);
		} else if ((s <= del) && (del < 2*s)) {
			remove_col(R.C, del - s);
			remove_col(R.D, del - s);
			remove_row(R.D, del - s);
			remove_row(R.gT, del - s);
		}
		remove_col(R.Z, del);
	}

	void remove_col(MatrixPartition<s> &A, int del) {
		int rows = A.rows();
		int cols = A.cols();

		if (del < cols) {
			A.block(0, del, rows, cols-del-1) = A.block(0, del+1, rows, cols-del-1);
		}
		A.conservativeResize(rows, cols-1);
	}

	void remove_col(BigMatrixPartition<s> &A, int del) {
		int rows = A.rows();
		int cols = A.cols();

		if (del < cols) {
			A.block(0, del, rows, cols-del-1) = A.block(0, del+1, rows, cols-del-1);
		}
		A.conservativeResize(rows, cols-1);
	}

	void remove_row(MatrixPartition<s> &A, int del) {
		int rows = A.rows();
		int cols = A.cols();

		if (del < rows) {
        	A.block(del, 0, rows-del-1, cols) = A.block(del+1, 0, rows-del-1, cols);
		}
		A.conservativeResize(rows-1, cols);
	}

	void remove_row(VectorPartition<s> &A, int del) {
		int rows = A.rows();
		int cols = A.cols();

		if (del < rows) {
        	A.block(del, 0, rows-del-1, cols) = A.block(del+1, 0, rows-del-1, cols);
		}
		A.conservativeResize(rows-1, cols);
	}

	/*! Concatenates the derivatives (parametrised beta and decision varibles f) into the RHS of (1b).
	
		\param B0 Parameterised derivatives at beginning of segement.
		\param f0 Decision variable derivatives at beginning of segment.
		\param BT Parameterised derivatives at end of segement.
		\param fT Decision variable derivatives at end of segment.
	*/
	VectorB<s> concatenate_derivatives(double B0, VectorPartition<s> f0, double BT, VectorPartition<s> fT) {
		VectorB<s> f;

		Scalar B0s; 
		B0s << B0;
		Scalar BTs; 
		BTs << BT;

		f << B0s, f0, BTs, fT;

		return f;
	}

	// Get and set functions.
	// NOTE: MinimisingPoly allocates and frees memory for P, d and Hi. It handles T and B using pointers, assuming whatever calls it allocates and frees the memory for T and B appropriately.

	/*! Sets time and derivative parameters. 

		\param T_data Array of segment times.   
		\param B_data Array of fixed derivative values.
	*/
	void set_vec_waypoint(double *T_data, double *B_data) {
		T = T_data;
		// B = B_data;
	}

	void set_map_waypoint(double *T_data, std::map<int, double> *B_data) {
		T = T_data;
		B = B_data;
	}

	/*! Returns pointer to array of vectors of polynomial coefficients of segments. 
	*/	
	VectorP<s> * get_P() {
        return P;
    }

	/*! Returns pointer to array of vectors of derivatives at the beginning and end of segments. 
	*/	
	VectorP<s> * get_d() {
		return d;
	}

	MatrixBlock<s> * get_H() {
		if (save_mem) {
			std::cout << "Warning, didn't save matrices. Returning pointer to nothing.";
		}
		return Hi;
	}	

	void refresh_block_matrices() {
		if (save_mem) {
			std::cout << "Warning, didn't save matrices. No memory prepared for matrices.";
		}
		gen_block_matrices(Hi);
	}
};

typedef MinimisingPoly<5> SnapMinimiser;
typedef MinimisingPoly<4> JerkMinimiser;
typedef MinimisingPoly<3> AcclMinimiser;

void saveTrajectory(double *t, std::map<int, double> *b, int k, std::string &str) {
    SnapMinimiser opt_traj(k, false, false);
	opt_traj.set_map_waypoint(t, b);

    VectorP<5> * P_opt = opt_traj.get_P();
    opt_traj.solve_minimisation();

    MatrixBlock<5> *Hi = opt_traj.get_H();
    double J = 0.0;
    for (int i = 0; i < k; i++) {
        J += pow(t[i+1] - t[i], -7.0) * (1.0 / 2) * P_opt[i].transpose() * Hi[i].Q * P_opt[i];
    }

    std::cout << std::fixed;
    std::ofstream myfile;
	
	double *bb;
	bb = new double [k+1];
	for (int i = 0; i < k; i++) {
		std::map<int, double>::iterator itr = b[i].begin();
		bb[i] = itr->second;
		if (i == k-1) {
			itr++;
			bb[k] = itr->second;
		}
	}

    myfile.open(str);
    myfile << k << '\n';
    for (int i = 0; i < k + 1; i++) {
        myfile << std::to_string(t[i]) << ' ';
    }
    myfile << '\n';
    for (int i = 0; i < k + 1; i++) {
        myfile << std::to_string(bb[i]) << ' ';
    }
    myfile << '\n';
    myfile << std::to_string(J) << '\n';
    for (int i = 0; i < k; i++) {		
        myfile << P_opt[i].transpose() << '\n';
    }
    myfile.close();
}

void save_jerk_trajectory(double *t, double *b, int k, std::string &str) {
    JerkMinimiser opt_traj(k, false, false);
    opt_traj.set_vec_waypoint(t, b);

    VectorP<4> * P_opt = opt_traj.get_P();
    opt_traj.solve_minimisation();

    MatrixBlock<4> *Hi = opt_traj.get_H();
    double J = 0.0;
    for (int i = 0; i < k; i++) {
        J += pow(t[i+1] - t[i], -7.0) * (1.0 / 2) * P_opt[i].transpose() * Hi[i].Q * P_opt[i];
    }

    std::cout << std::fixed;
    std::ofstream myfile;

    myfile.open(str);
    myfile << k << '\n';
    for (int i = 0; i < k + 1; i++) {
        myfile << std::to_string(t[i]) << ' ';
    }
    myfile << '\n';
    for (int i = 0; i < k + 1; i++) {
        myfile << std::to_string(b[i]) << ' ';
    }
    myfile << '\n';
    myfile << std::to_string(J) << '\n';
    for (int i = 0; i < k; i++) {		
        myfile << P_opt[i].transpose() << '\n';
    }
    myfile.close();
}

void generateTrajectory(int k, double* t, double *bb) {
    std::random_device rd;
    std::mt19937 gen(rd());
    // std::mt19937 gen(3);

    t[0] = 0.0;
    bb[0] = (double(gen() % 2000) / 1000.0) - 1;
    for (int i = 1; i < k + 1; i ++) {
        t[i] = t[i-1] + (double(gen() % 1000) / 1000.0) + 0.1;
        bb[i] = (double(gen() % 2000) / 1000.0) - 1 + bb[i-1];
    }
}