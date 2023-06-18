#include <algorithm>

#include "MathModel.h"



void MathModel::DriftModel::SetInitialConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	const std::valarray<double>& theta,
	double dz)
{
	switch (_task_type)
	{
	case MathModel::TaskType::BubblesRising:
		SetBubblesRisingInitialConditions(alpha_g, p, v_m, v_g, v_l, theta, dz);
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
	std::valarray<double>C_0(1.0, n_points_cell_velocities);

	double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
	double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
	double dz_dimensionless = dz / _L; // ������������ ��� �� ������������

	for (size_t i = 0; i < n_points_cell_velocities; ++i)
	{
		double R = (d[i] + d[i + 1]) / 4; // ������ �����

		// ������������� ����������
		const double pi = 3.14;
		double r_0 = 1 * R; // ���������� ����� ����������
		double x_0 = r_0 / R; // ���� ������� �����, �� ������� ����� ���
		double mu = GetMixtureViscosity(_alpha_g_0); // ����������� �������� �����
		const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
		double r_b = r_b0; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
		double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
		double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
		double St = mb_0 * _U  / (6 * pi * mu_0 * R * r_b0); // ����� C�����
		double Fr = _U / sqrt(g * R); // ����� �����
		double Re = rho_l_0 * _U * R / mu_0;
		double eta = rho_g_0 / rho_l_0; // ��������� ��������� ���� � ��������� ��������
		double rho_l = GetLiquidDensity(p[i]) / rho_l_0; //������������ ������� ��������� ��������
		double rho_g = GetGasDensity(p[i], dz * i) / rho_g_0; // ������������ ������� ��������� ����
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

	return C_0 * _U;
}

std::valarray<double> MathModel::DriftModel::CalculateV_d(const std::valarray<double>& d, const std::valarray<double>& alpha_g, std::valarray<double> p, size_t n_points_cell_velocities)
{
	std::valarray<double>v_d(n_points_cell_velocities);

	for (size_t i = 0; i < n_points_cell_velocities; ++i)
	{
		const double pi = 3.14;
		double R = (d[i] + d[i + 1]) / 4; // ������ �����
		double mu = GetMixtureViscosity(_alpha_g_0); // ����������� �������� �����
		const double r_b0 = 1E-3; // ����������� �������� ��� ������� ��������
		double r_b = 1; // ��� ������ ����������� ����� ���������. ������ ���� ��������� �������
		double rho_l_0 = GetCharacteristicLiquidDensity(); // ����������� ��������� ��������
		double rho_g_0 = GetCharacteristicGasDensity(); // ����������� ��������� ����
		double mb_0 = 4 / 3 * pi * pow(r_b0, 3) * rho_g_0; // ����� ��������
		double mu_0 = GetLiquidViscosity(); // �������� ������ ��������
		double St = mb_0 * _U  / (6 * pi * mu_0 * R * r_b0); // ����� C�����
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

double MathModel::DriftModel::GasSteadyFlowAnalyticsVelocity(double p, double z, double F)
{
	double R = 8.314; // ������� ����������
	double T = 293;   // ����������� � ���������
	double C1 = -0.03567;
	double v = C1 * R * T / (p * F);
	return v;
}

double MathModel::DriftModel::GasSteadyFlowAnalyticsDensity(double v, double p)
{
	double R = 8.314; // ������� ����������
	double T = 293;   // ����������� � ���������
	double C2 = 28076.974;

	double must_be_zero = 0.5 * v * v + R * T * log(p) - C2;

	
	return must_be_zero;
}

double MathModel::DriftModel::CalculateHydrostaticPressure(double rho, double h)
{
	return rho * g * h;
}


double MathModel::DriftModel::GetMixtureDensity(double alpha_g, double alpha_l, double p, double z)
{
	return alpha_g * GetGasDensity(p, z) + alpha_l * GetLiquidDensity(p);
}

double MathModel::DriftModel::GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l)
{
	return alpha_g * v_g + alpha_l * v_l;
}

double MathModel::DriftModel::GetMixtureViscosity(double alpha_g)
{
	return 1 / (1 - alpha_g);
}

double MathModel::DriftModel::GetGasDensity(double p, double z)
{
	
	// return  GetCharacteristicGasDensity(); // ��� ������ ����������� ������������ ���������
	double R = 8.314; // ������� ����������
	double T = 293;   // ����������� � ���������
	return p / (R * T);
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
	double water_viscosity = 1; // ������� �������� ������ �������� (����� � ������������ ����)
	return water_viscosity;
}

double MathModel::DriftModel::GetFrictionCoefficient(double alpha_g,double alpha_l, double v_m, double p, double d, double eps)
{
	double rho_m = GetMixtureDensity(alpha_g, alpha_l, p, 0); // �� ������ ����������� ���������
	double mu_m = GetMixtureViscosity(alpha_g);
	v_m = abs(v_m);
	double Re = (v_m * d * rho_m) / (mu_m * 1E-3);
	double f;

	if (Re > 2320)
	{
		f = 9.09 * pow(rho_m * v_m * d / (eps * rho_m * v_m + 68 * mu_m), 0.25); // ������� �������� (155 �������� � "����������� ������������ ��������� ����")
	}
	else 
	{

		f = Re != 0 ? 64 / Re : 0;
	}


	return f;
}

void MathModel::DriftModel::SetBubblesRisingInitialConditions(
	std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	const std::valarray<double>& theta,
	double dz)
{
	for (size_t i = 0; i < alpha_g.size(); ++i)
	{
		alpha_g[i] = 0;
	}

	for (size_t i = 0; i < p.size(); ++i)
	{
		p[i] = atm + CalculateHydrostaticPressure(GetCharacteristicLiquidDensity(), (i) * dz) * cos(theta[i]);
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
	double gas_flow = GetBubblesRisingGasFlow();
	double liquid_flow = GetBubblesRisingLiquidFlow();
	double S = well.GetBottomCrossSectionArea();

	// ����� ��������
	size_t index_wb_velocity = v_g.size() - 1;
	size_t index_wb_property = alpha_g.size() - 1;

	// ����� ��������
	size_t index_wt = 0;

	// ������� �� �������
	p[index_wt] = atm;

	// ������� �� �������� ����
	alpha_g[index_wb_property] = gas_flow / (gas_flow + liquid_flow);
	
	
	// ������� �� ��������
	v_g[index_wb_velocity] = - gas_flow / S;
	v_l[index_wb_velocity] = - liquid_flow / S;
	v_m[index_wb_velocity] = GetMixtureVelocity(
		(alpha_g[index_wb_property] + alpha_g[index_wb_property - 1]) / 2,
	 1 - (alpha_g[index_wb_property] + alpha_g[index_wb_property - 1]) / 2,
			v_g[index_wb_velocity],
			v_l[index_wb_velocity]
		);
}

void MathModel::DriftModel::SetDebugInitialConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, std::valarray<double>& v_g, std::valarray<double>& v_l, double dz)
{
	for (size_t i = 0; i < alpha_g.size(); ++i)
	{
		alpha_g[i] = 0;
	}

	for (size_t i = 0; i < p.size(); ++i)
	{
		p[i] = (CalculateHydrostaticPressure(GetCharacteristicLiquidDensity(), (i + 1) * dz));
	}

	for (size_t i = 0; i < v_m.size(); ++i)
	{
		v_l[i] = 0;
		v_g[i] = 0;
		v_m[i] = 0;
	}
}



double MathModel::DriftModel::GetBubblesRisingLiquidFlow()
{
	double flow_value = 0.1 / 600; // 100 ������ / ������
	return flow_value * 10;
}

double MathModel::DriftModel::GetBubblesRisingGasFlow()
{
	double flow_value = GetBubblesRisingLiquidFlow() / 19;
	return flow_value;
}

void MathModel::DriftModel::SetBubblesRisingCharacteristicVelocity(Well well)
{
	double gas_velocity = GetBubblesRisingGasFlow() / well.GetBottomCrossSectionArea();
	double liquid_velocity = GetBubblesRisingLiquidFlow() / well.GetBottomCrossSectionArea();
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

void MathModel::DriftModel::SetDebugBoundaryConditions(std::valarray<double>& alpha_g,
	std::valarray<double>& p,
	std::valarray<double>& v_m,
	std::valarray<double>& v_g,
	std::valarray<double>& v_l,
	const Well well,
	double dt)
{

	size_t index_velocity = 0;
	size_t index_property = 1;

	// ������� �� �������� ����
	// alpha_g[index_property] = 0.05;

	// ������� �� ��������
	// v_g[index_velocity] = 0.2;
	// v_l[index_velocity] = 0;

	// �������� �������� �� ������ ��� ��������
	double alpha_g_mid = (alpha_g[index_property] + alpha_g[index_property - 1]) / 2;
	double alpha_l_mid = 1 - alpha_g_mid;


	// v_m[index_velocity] = GetMixtureVelocity(alpha_g_mid, alpha_l_mid, v_g[index_velocity], v_l[index_velocity]);

}