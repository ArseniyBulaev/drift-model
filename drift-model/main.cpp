#include <iostream>
#include <string>
#include <fstream>
#include <valarray>

#include "Well\Well.h"
#include "Well\WellSegment.h"
#include "Solver\DriftModelSolver.h"

using namespace std;


void write_to_file(const std::valarray<double> & v, const double & dz, const string file_name);

int main() {

	double dz = 10;
	double dt = 0.00005;
	Well well = Well({ WellSegment(3000,0,0,0.1)});
	MathModel::TaskType task_type = MathModel::TaskType::BubblesRising;
	DriftModelSolver solver(dz, dt, well, task_type);

	solver.Solve();

	const double atm = 98'066.5;

	write_to_file(solver.Get_P(), solver.GetDz(), "p.txt");
	write_to_file(solver.GetV_m(), solver.GetDz(), "v_m.txt");
	write_to_file(solver.GetV_g(), solver.GetDz(), "v_g.txt");
	write_to_file(solver.GetAlpha_g(), solver.GetDz(), "alpha_g.txt");
	

	return 0;
}


void write_to_file(const std::valarray<double> & v, const double & dz, const string file_name) {

	std::string files_dirrectory = "Results";
	std::ofstream myfile(files_dirrectory + "/" + file_name);
	for (size_t i = 0; i < v.size(); ++i)
	{
		myfile << dz * (i + 1) << "\t" << v[i] << std::endl;
	}
	myfile.close();
}