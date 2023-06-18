#include <iostream>
#include <string>
#include <fstream>
#include <valarray>

#include "Well\Well.h"
#include "Well\WellSegment.h"
#include "Solver\DriftModelSolver.h"

using namespace std;


int main() {


	double dz = 1;
	double dt = 0.001;
	double T = 500; // Время расчёта в секундах
	Well well = Well({ WellSegment(10,0,0,0.1)});
	MathModel::TaskType task_type = MathModel::TaskType::BubblesRising;


	// DriftModelSolver solver(T, dz, dt, well, task_type);
	// DriftModelSolver solver(T, dz / 10, dt, well, task_type);
	// DriftModelSolver solver(T, dz / 100, dt, well, task_type);
	DriftModelSolver solver(T, dz / 1000, dt, well, task_type);


	solver.Solve();



	return 0;
}


