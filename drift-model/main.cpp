#include <iostream>
#include <string>
#include <fstream>
#include <valarray>

#include "Well\Well.h"
#include "Well\WellSegment.h"
#include "Solver\DriftModelSolver.h"

using namespace std;


int main() {

	double dz = 0.001;
	double dt = 0.000001;
	Well well = Well({ WellSegment(10,0,0,0.1)});
	MathModel::TaskType task_type = MathModel::TaskType::BubblesRising;
	DriftModelSolver solver(dz, dt, well, task_type);

	solver.Solve();
	

	return 0;
}


