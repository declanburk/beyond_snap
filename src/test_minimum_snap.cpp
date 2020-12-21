#include <iostream>
#include <cmath>
#include <fstream>
#include <map>

#include <minimising_poly.h>
#include <timer.h>

#define N 6

using std::fixed;

void oneLap(int K, bool save) {
	// Object to generate minimum snap trajectories.
	bool save_mem = true;
	bool debug = false;
	SnapMinimiser test_min(K, save_mem, debug);

	// generate random walk to interpolate. This for() loop is responsible for T and B memory. SnapMinimiser takes care of P and d.
	double *T, *B;
	T = new double[K+1];
	B = new double[K+1];
	generateTrajectory(K, T, B);

	std::map<int, double> *Bmap;
	Bmap = new std::map<int, double> [K];
	for (int i = 0; i < K; i++) {
		Bmap[i].insert(std::pair<int, double>(0, B[i]));
		if (i == 0) {
			Bmap[i].insert(std::pair<int, double>(1, 0.0));
		}
		Bmap[i].insert(std::pair<int, double>(5, B[i+1]));
		if (i == K-1) {
			Bmap[i].insert(std::pair<int, double>(6, 0.0));
		}
	}

	// test_min.set_vec_waypoint(T, B);
	test_min.set_map_waypoint(T, Bmap);
	test_min.solve_minimisation();

	if (save) {
		std::string name = "min_snap_coefficients.csv";
		saveTrajectory(T, Bmap, K, name);
	}

	delete[] T;
	delete[] B;
	delete[] Bmap;
}

int main() {
	// run (N) simulations and repeat each one (sim_repetitions) times.
	int sim_repetitions = 1;
	int sim_waypoints [N] = {10, 50, 100, 500, 1000, 5000};
	double sim_timing [N] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double total_sim_timing [N] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	// timing variables
	clock_t start, end;
    double cpu_time_used, total_cpu_time_used;
	for (int j = 0; j < N; j++) {
		int K = sim_waypoints[j]; // number of segements

		LapTimer tim;
		total_cpu_time_used = 0;
		for (int i = 0; i < sim_repetitions; i++) {
			tim.start();
			oneLap(K, false);			
			total_cpu_time_used += tim.stop();
		}

		total_sim_timing[j] = total_cpu_time_used;
		sim_timing[j] = tim.avg;

		std::cout << "Calculated trajectory of " << sim_waypoints[j] << " segments in " << total_sim_timing[j] << " s. Average time is " << sim_timing[j] << " s.\n";
	}

	oneLap(sim_waypoints[5], true);

	return 0;
}