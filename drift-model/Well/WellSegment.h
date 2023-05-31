#pragma once
class WellSegment
{

public:
	const double length = 0.0; // Длина
	const double tilt_angle = 0.0; // Угол наклона 
	const double relative_roughness = 0.0; // Относительная шероховатость
	const double diameter = 0.0; // Диаметр

	WellSegment(double l, double theta, double eps, double d);
	
};

