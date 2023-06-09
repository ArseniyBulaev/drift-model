#include <cmath>

#include "MathModel.h"



void MathModel::DriftModel::SetInitialConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	double dz)
{
	switch (_task_type)
	{
	case MathModel::TaskType::BubblesRising:
		SetBubblesRisingInitialConditions(alpha_g, p, v_m, v_g, v_l, dz);
		break;
	case MathModel::TaskType::Debug:
		SetDebugInitialConditions(alpha_g, p, v_m, v_g, v_l, dz);
		break;
	}
}

void MathModel::DriftModel::SetBoundaryConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	const Well well,
	double dt)
{
	switch (_task_type)
	{
	case MathModel::TaskType::BubblesRising:
		SetBubblesRisingBoundaryConditions(alpha_g, p, v_m , v_g, v_l, well, dt);
		break;
	case MathModel::TaskType::Debug:
		SetDebugBoundaryConditions(alpha_g, p, v_m, v_g, v_l, well, dt);
		break;
	}
}

std::valarray<double> MathModel::DriftModel::CalculateC_0(const std::valarray<double>& d, const std::valarray<double>& alpha_g, const std::valarray<double>& p, size_t n_points_cell_velocities, double dz)
{
	std::valarray<double>C_0(n_points_cell_velocities);

	double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
	double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
	double dz_dimensionless = dz / _L; // ������������ ��� �� ������������

	for (size_t i = 0; i < n_points_cell_velocities; ++i)
	{
		double R = (d[i] + d[i + 1]) / 4; // ������ �����
		_alpha_g_0 = (alpha_g[i] + alpha_g[i + 1]) / 2;

		// ������������� ����������
		const double pi = 3.14;
		double r_0 = 1 * R; // ���������� ����� ����������
		double x_0 = r_0 / R; // ���� ������� �����, �� ������� ����� ���
		double mu = GetMixtureViscosity(_alpha_g_0); // ����������� �������� �����
		const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
		double r_b = r_b0; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
		double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
		double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
		double St = mb_0 * _U * pow(r_b, 3) / (6 * pi * mu_0 * R * r_b0); // ����� C�����
		double Fr = _U / sqrt(g * R); // ����� �����
		double Re = rho_l_0 * _U * R / mu_0;
		double eta = rho_g_0 / rho_l_0; // ��������� ��������� ���� � ��������� ��������
		double rho_l = GetLiquidDensity(p[i]) / rho_l_0; //������������ ������� ��������� ��������
		double rho_g = GetGasDensity(p[i]) / rho_g_0; // ������������ ������� ��������� ����
		double eps = R / _L; // ����� �������� (����������� �������� ������)

		double dp_dz = eps * Re * (p[i + 1] - p[i]) / dz_dimensionless; // ������������ �������� ��������
		

		// ���������� � �������
		double m_0 = 1 / (r_b * mu);
		double gamma = St / (eta * pow(Fr, 2));
		double a = (pow(R, 2) / 4) * (-dp_dz + (Re / pow(Fr, 2)) * rho_l);
		double b = (pow(R, 2) / 4) * (Re / pow(Fr, 2)) * (rho_l - eta * rho_g) * _alpha_g_0;

		// �������� ��������� �������
		C_0[i] = (2 * a * (1 - pow(x_0, 2)) + m_0 * (a - b) * pow(x_0, 2) - 2 * gamma * rho_l * m_0 * _alpha_g_0) / (a * (1 - pow(x_0, 4)) + m_0 * (a - b) * pow(x_0, 4) - 2 * gamma * rho_l * m_0 * _alpha_g_0 * pow(x_0, 2));
	}

	return C_0;
}

std::valarray<double> MathModel::DriftModel::CalculateV_d(const std::valarray<double>& d, const std::valarray<double>& alpha_g, std::valarray<double> p, size_t n_points_cell_velocities)
{
	std::valarray<double>v_d(n_points_cell_velocities);

	for (size_t i = 0; i < n_points_cell_velocities; ++i)
	{
		_alpha_g_0 = (alpha_g[i] + alpha_g[i + 1]) / 2;
		const double pi = 3.14;
		double R = (d[i] + d[i + 1]) / 4; // ������ �����
		double mu = GetMixtureViscosity(_alpha_g_0); // ����������� �������� �����
		const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
		double r_b = r_b0; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
		double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
		double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
		double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
		double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
		double St = mb_0 * _U * pow(r_b, 3) / (6 * pi * mu_0 * R * r_b0); // ����� C�����
		double Fr = _U / sqrt(g * R); // ����� �����
		double eta = rho_g_0 / rho_l_0; // ��������� ��������� ���� � ��������� ��������
		double rho_l = GetLiquidDensity(p[i]) / rho_l_0; // ������� ��������� ��������

		// ���������� � �������
		double m_0 = 1 / (r_b * mu);
		double gamma = St / (eta * pow(Fr, 2));

		// �������� ������
		v_d[i] = -gamma * rho_l * m_0 * (1 - _alpha_g_0);
	}

	return v_d * _U;
}

double MathModel::DriftModel::CalculateHydrostaticPressure(double rho, double h)
{
	return rho * g * h;
}



double MathModel::DriftModel::CalculateGasProfileParameter_TEST(double alpha_g)
{
	const double C_0_0 = 1.2;
	const double b = 0.6;

	double C_0 = C_0_0 * (alpha_g < b ? 1 : 1.0 / (1 + (C_0_0 - 1) * pow((alpha_g - b) / (1 - b), 2)));

	return C_0;
}

double MathModel::DriftModel::CalculateDriftVelocity_TEST(double alpha_g)
{
	const double c_l = 1500; // �������� ����� � ����
	const double M_d = 1.5 * 1E-4;
	const double b_1 = 0.9;

	double v_d = c_l * M_d * (alpha_g < b_1 ? 1 : 1 - pow((alpha_g - b_1) / (1 - b_1), 2));

	return v_d;
}

double MathModel::DriftModel::GetMixtureDensity(double alpha_g, double alpha_l, double p)
{
	return alpha_g * GetGasDensity(p) + alpha_l * GetLiquidDensity(p);
}

double MathModel::DriftModel::GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l)
{
	return alpha_g * v_g + alpha_l * v_l;
}

double MathModel::DriftModel::GetMixtureViscosity(double alpha_g)
{
	return 1 / (1 - alpha_g);
}

double MathModel::DriftModel::GetGasDensity(double p)
{
	return  GetCharacteristicGasDensity(); // ��� ������ ����������� ������������ ���������
}

double MathModel::DriftModel::GetLiquidDensity(double p)
{
	
	return GetCharacteristicLiquidDensity(); // ��� ������ ����������� ������������ ���������
	
}

double MathModel::DriftModel::GetCharacteristicGasDensity()
{
	double density = 1.29; // ��������� �������
	return density;
}

double MathModel::DriftModel::GetCharacteristicLiquidDensity()
{
	double density = 997; // ��������� ����
	return density;
}



double MathModel::DriftModel::GetLiquidViscosity()
{
	double water_viscosity = 1; // ������� �������� ������ ��������
	return water_viscosity;
}

double MathModel::DriftModel::GetFrictionCoefficient(double alpha_g,double alpha_l, double v_m, double p, double d, double eps)
{
	double rho_m = GetMixtureDensity(alpha_g, alpha_l, p);
	double mu_m = GetMixtureViscosity(alpha_g);


	double f = v_m > 0 ? 9.09 * pow(rho_m * v_m * d / (eps * rho_m * v_m + 68 * mu_m), 0.25) : 1; // ������� �������� (155 �������� � "����������� ������������ ��������� ����") 

	return 0;
}

void MathModel::DriftModel::SetBubblesRisingInitialConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	double dz)
{
	for (size_t i = 0; i < alpha_g.size(); ++i)
	{
		alpha_g[i] = 0;
	}

	for (size_t i = 0; i < p.size(); ++i)
	{
		p[i] = CalculateHydrostaticPressure(GetCharacteristicLiquidDensity(), (i + 1) * dz);
	}
	
	for (size_t i = 0; i < v_m.size(); ++i)
	{
		v_l[i] = 0;
		v_g[i] = 0;
		v_m[i] = 0;
	}

}

void MathModel::DriftModel::SetBubblesRisingBoundaryConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	const Well well,
	double dt)
{
	double gas_flow = GetBubblesRisingGasFlow(dt);
	double liquid_flow = GetBubblesRisingLiquidFlow(dt);
	double S = well.GetBottomCrossSectionArea();

	// ��������� ������� �������� �� �����
	size_t index_wb_velocity = v_g.size() - 1;
	size_t index_wb_property = alpha_g.size() - 1;

	// ������� �� �������� ����
	alpha_g[index_wb_property] = gas_flow / (gas_flow + liquid_flow);
	
	// ������� �� ��������
	v_g[index_wb_velocity] = - gas_flow / S;
	v_l[index_wb_velocity] = - liquid_flow / S;

	// �������� �������� �� ������ ��� ��������
	double alpha_g_mid = (alpha_g[index_wb_property] + alpha_g[index_wb_property - 1]) / 2;
	double alpha_l_mid = 1 - alpha_g_mid;


	v_m[index_wb_velocity] = GetMixtureVelocity(alpha_g_mid, alpha_l_mid, v_g[index_wb_velocity], v_l[index_wb_velocity]);
}

void MathModel::DriftModel::SetDebugInitialConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, std::valarray<double>& v_g, std::valarray<double>& v_l, double dz)
{
	for (size_t i = 0; i < alpha_g.size(); ++i)
	{
		alpha_g[i] = 0;
	}

	for (size_t i = 0; i < p.size(); ++i)
	{
		p[i] = CalculateHydrostaticPressure(GetCharacteristicLiquidDensity(), (i + 1) * dz);
	}

	for (size_t i = 0; i < v_m.size(); ++i)
	{
		v_l[i] = 0;
		v_g[i] = 0;
		v_m[i] = 0;
	}
}

void MathModel::DriftModel::SetDebugBoundaryConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, std::valarray<double>& v_g, std::valarray<double>& v_l, const Well well, double dt)
{
	
	size_t index_velocity = 0;
	size_t index_property = 1;

	// ������� �� �������� ����
	alpha_g[index_property] = 0.05;

	// ������� �� ��������
	v_g[index_velocity] = 0.2;
	v_l[index_velocity] = 0;

	// �������� �������� �� ������ ��� ��������
	double alpha_g_mid = (alpha_g[index_property] + alpha_g[index_property - 1]) / 2;
	double alpha_l_mid = 1 - alpha_g_mid;


	v_m[index_velocity] = GetMixtureVelocity(alpha_g_mid, alpha_l_mid, v_g[index_velocity], v_l[index_velocity]);

}

double MathModel::DriftModel::GetBubblesRisingLiquidFlow(double dt)
{
	double flow_value = 1.0 / 3600; // 1 �^3 / ���
	return flow_value * dt;
}

double MathModel::DriftModel::GetBubblesRisingGasFlow(double dt)
{
	double flow_value = GetBubblesRisingLiquidFlow(dt) / 5;
	return flow_value;
}

void MathModel::DriftModel::SetBubblesRisingCharacteristicVelocity(Well well)
{
	double gas_velocity = GetBubblesRisingGasFlow(1) / well.GetBottomCrossSectionArea();
	double liquid_velocity = GetBubblesRisingLiquidFlow(1) / well.GetBottomCrossSectionArea();
	_U =  std::max(gas_velocity, liquid_velocity);
}


void MathModel::DriftModel::SetDebugCharacteristicVelocity(Well well)
{
	_U =  0.2;
}

void MathModel::DriftModel::SetDebugCharacteristicGasVolumeFraction()
{
	_alpha_g_0 = 0;
}
