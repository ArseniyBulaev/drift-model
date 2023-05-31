#include <stdexcept>
#include <cmath>

#include "DriftModel.h"


double DriftModel::CalculateHydrostaticPressure(double rho, double h)
{
	return rho * g * h;
}

double DriftModel::CalculateGasProfileParameter(double alpha_g_0, double U, double d,  double p, double dp_dz)
{
	const double pi = 3.14;
	const double g = 9.81; // ������ ���� ����������, ��� ��� ����� ���� ����������� ����� �����. � ������ ���� �������� ���� g, �� ����� ������ ������� �� ����

	// TO DO: ��������� ����� ���������� � ���������� ���� (������������)
	// ������������� ����������
	double R = d / 2; // ������ �����
	double r_0 = 0; // ��������� ���
	double x_0 = r_0 / R; // ���� �����, �� ������� ����� ���
	double mu = GetMixtureViscosity(alpha_g_0); // ����������� �������� �����
	const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
	double r_b = r_b0; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
	double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
	double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
	double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
	double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
	double St = mb_0 * U * pow(r_b, 3) / (6 * pi * mu_0 * R * r_b0); // ����� C�����
	double Fr = U / sqrt(g * R); // ����� �����
	double Re = rho_l_0 * U * R / GetLiquidViscosity();
	double eta = rho_g_0 / rho_l_0; // ��������� ��������� ���� � ��������� ��������
	double rho_l = GetLiquidDensity(p); // ������� ��������� ��������
	double rho_g = GetGasDensity(p); // ������� ��������� ����

	// ���������� � �������
	double m_0 = 1 / (r_b * mu);
	double gamma = St / (eta * pow(Fr, 2));
	double a = (pow(R, 2) / 4) *(-dp_dz + (Re / pow(Fr, 2)) * rho_l);
	double b = (pow(R, 2) / 4) * (Re / pow(Fr, 2)) * (rho_l - eta * rho_g) * alpha_g_0;

	// �������� ��������� �������
	double C_0 = (2*a*(1 - pow(x_0, 2)) + m_0 * (a - b) * pow(x_0, 2) - 2 * gamma * rho_l * m_0 * alpha_g_0) / (a * (1 - pow(x_0, 4)) + m_0 * (a - b) * pow(x_0, 4) - 2 * gamma * rho_l * m_0 * alpha_g_0 * pow(x_0, 2));

	return 2;
	return C_0;
}

double DriftModel::CalculateDriftVelocity(double alpha_g_0, double U, double d, double p)
{
	const double pi = 3.14;
	const double g = 9.81; // ������ ���� ����������, ��� ��� ����� ���� ����������� ����� �����. � ������ ���� �������� ���� g, �� ����� ������ ������� �� ����

	// TO DO: ��������� ����� ���������� � ���������� ���� (������������)
	// ������������� ����������
	double R = d / 2; // ������ �����
	double mu = GetMixtureViscosity(alpha_g_0); // ����������� �������� �����
	const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
	double r_b = r_b0; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
	double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
	double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
	double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
	double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
	double St = mb_0 * U * pow(r_b, 3) / (6 * pi * mu_0 * R * r_b0); // ����� C�����
	double Fr = U / sqrt(g * R); // ����� �����
	double eta = rho_g_0 / rho_l_0; // ��������� ��������� ���� � ��������� ��������
	double rho_l = GetLiquidDensity(p); // ������� ��������� ��������

	// ���������� � �������
	double m_0 = 1 / (r_b * mu);
	double gamma = St / (eta * pow(Fr, 2));

	// �������� ������
	double v_d = -gamma * rho_l * m_0 * (1 - alpha_g_0);

	return 0;
	return v_d;
}

double DriftModel::CalculateGasProfileParameter_TEST(double alpha_g)
{
	const double C_0_0 = 1.2;
	const double b = 0.6;

	double C_0 = C_0_0 * (alpha_g < b ? 1 : 1.0 / (1 + (C_0_0 - 1) * pow((alpha_g - b) / (1 - b), 2)));

	return C_0;
}

double DriftModel::CalculateDriftVelocity_TEST(double alpha_g)
{
	const double c_l = 1500; // �������� ����� � ����
	const double M_d = 1.5 * 1E-4;
	const double b_1 = 0.9;

	double v_d = c_l * M_d * (alpha_g < b_1 ? 1 : 1 - pow((alpha_g - b_1) / (1 - b_1), 2));

	return v_d;
}

double DriftModel::GetMixtureDensity(double alpha_g, double alpha_l, double p)
{
	return alpha_g * GetGasDensity(p) + alpha_l * GetLiquidDensity(p);
}

double DriftModel::GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l)
{
	return alpha_g * v_g + alpha_l * v_l;
}

double DriftModel::GetMixtureViscosity(double alpha_g)
{
	return 1 / (1 - alpha_g);
}

double DriftModel::GetGasDensity(double p)
{
	return  GetCharacteristicGasDensity(); // ��� ������ ����������� ������������ ���������
}

double DriftModel::GetLiquidDensity(double p)
{
	
	return GetCharacteristicLiquidDensity(); // ��� ������ ����������� ������������ ���������
	
}

double DriftModel::GetCharacteristicGasDensity()
{
	double density = 1.29; // ��������� �������
	return density;
}

double DriftModel::GetCharacteristicLiquidDensity()
{
	double density = 1000.0; // ��������� �������������� ��������
	return density;
}

double DriftModel::GetGasViscosity()
{
	double air_viscosity = 18.35E-6;
	return air_viscosity;
}

double DriftModel::GetLiquidViscosity()
{
	double water_viscosity = 8.9E-4; // �������� ���� ��� ����������� 25 �������� �� �������
	return water_viscosity;
}


double DriftModel::GetFrictionCoefficient(double alpha_g,double alpha_l, double v_m, double p, double d, double eps)
{
	double rho_m = GetMixtureDensity(alpha_g, alpha_l, p);
	double mu_m = GetMixtureViscosity(alpha_g);


	double f = v_m > 0 ? 9.09 * pow(rho_m * v_m * d / (eps * rho_m * v_m + 68 * mu_m), 0.25) : 1; // ������� �������� (155 �������� � "����������� ������������ ��������� ����") 

	return 1;
}


double DriftModel::GetGSGasVolumeFractionBoundaryCondition()
{
	return 0.0;
}

double DriftModel::GetGSPressureBoundaryCondition()
{
	return 0.0;
}

double DriftModel::GetGSLiquidVelocityBoundaryCondition()
{
	return 0.0;
}

double DriftModel::GetGSGasVelocityBoundaryCondition()
{
	return 0.0;
}

double DriftModel::GetGSMixtureVelocityBoundaryCondition()
{
	return 0.0;
}

double DriftModel::GetGSGasVolumeFractionInitialCondition(double L, double z)
{
	return z > L / 2 ? 1 : 0.0;
}

double DriftModel::GetPRCharacteristicVelocity(double p_wf, double well_bottom_cross_section_area)
{
	double p_r = 3.5E+7; // ��������� ��������
	double PI = 4.4 * 1E-11;
	double q = PI * (p_r - p_wf);
	return q / well_bottom_cross_section_area;
}

double DriftModel::GetPRCharacteristicStopTime(double p_wf, double well_bottom_cross_section_area, double well_length)
{
	double U = GetPRCharacteristicVelocity(p_wf, well_bottom_cross_section_area);
	return well_length / U;
}

double DriftModel::GetPRLiquidVelocityBoundaryCondition(double p_wf, double well_bottom_cross_section_area)
{
	return GetPRCharacteristicVelocity(p_wf, well_bottom_cross_section_area);
	
}


double DriftModel::GetPRLiquidVolumeFractionBoundaryCondition(double v_m, double v_l) 
{
	return v_m / v_l;
}

double DriftModel::GetPRGasVolumeFractionInitialCondition()
{
	return 0.0;
}
