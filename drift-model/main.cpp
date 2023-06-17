#include <iostream>
#include <string>
#include <fstream>
#include <valarray>

#include "Well\Well.h"
#include "Well\WellSegment.h"
#include "Solver\DriftModelSolver.h"

using namespace std;


int main() {


	double dz = 10;
	double dt = 1;
	double T = 566038; // Время расчёта в секундах
	Well well = Well({ WellSegment(3000,0,0,0.1)});
	MathModel::TaskType task_type = MathModel::TaskType::BubblesRising;
	DriftModelSolver solver(T, dz, dt, well, task_type);

	solver.Solve();
	

	return 0;
}


