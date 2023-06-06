#include <iostream>
#include <string>
#include <fstream>
#include <valarray>

#include "Well\Well.h"
#include "Well\WellSegment.h"
#include "Solver\DriftModelSolver.h"

using namespace std;


int main() {

	double dz = 0.01;
	double dt = 0.001;
	Well well = Well({ WellSegment(1,0,0,0.1)});
	MathModel::TaskType task_type = MathModel::TaskType::BubblesRising;
	DriftModelSolver solver(dz, dt, well, task_type);

	solver.Solve();
	

	return 0;
}


