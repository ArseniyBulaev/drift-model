#pragma once
class WellSegment
{

public:
	const double length = 0.0; // �����
	const double tilt_angle = 0.0; // ���� ������� 
	const double relative_roughness = 0.0; // ������������� �������������
	const double diameter = 0.0; // �������

	WellSegment(double l, double theta, double eps, double d);
	
};

