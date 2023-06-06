#include "Writer.h"


void Writer::WriteToFile(const std::valarray<double>& v, const double& dz, const std::string file_name)
{
	
	std::ofstream myfile(files_dirrectory + "/" + file_name);
	for (size_t i = 0; i < v.size(); ++i)
	{
		myfile << dz * i << "\t" << v[i] << std::endl;
	}
	myfile.close();
}

