#include <iostream>
#include <fstream>
#include <string>
#include <valarray>


#pragma once
class Writer
{
public:
	void WriteToFile(const std::valarray<double>& v, const double& dz, const std::string file_name);
private:
	const std::string files_dirrectory = "Results";
};

