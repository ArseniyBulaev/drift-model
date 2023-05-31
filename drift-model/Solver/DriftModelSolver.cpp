#include "DriftModelSolver.h"

//DEBUG
#include <iostream>



DriftModelSolver::DriftModelSolver(double dz, double dt,  const Well & well, MathModel::TaskType task_type): _dz(dz), _dt(dt), _well(well), _drift_model(task_type)
{
	// ���������� ����� ����� ��� �������� ���� ���������� ����� ��������� �� ��������� ���� dz
	_n_points_cell_properties = CalculateN();
	_n_points_cell_velocities = _n_points_cell_properties - 1;

	// ������������� ��������
	// ��������� �����
	_theta = std::valarray<double>(_n_points_cell_properties); // ���� ������� �����
	_d = std::valarray<double>(_n_points_cell_properties); // ������� �����
	_eps = std::valarray<double>(_n_points_cell_properties); // ������������� �������������
	// ��������� ���������
	_v_m = std::valarray<double>(_n_points_cell_velocities); // ������� �����
	_v_g = std::valarray<double>(_n_points_cell_velocities); // ������� ����
	_v_l = std::valarray<double>(_n_points_cell_velocities); // ������� ��������
	_alpha_g = std::valarray<double>(_n_points_cell_properties); // �������� ���� ����
	_alpha_l = std::valarray<double>(_n_points_cell_properties); // �������� ���� ��������
	_p = std::valarray<double>(_n_points_cell_properties); // �������� � ���������� �����


	// ������������� ���������� �����
	InitializeGeometryParameters();

	// ������������� ��������� �������
	initializeCellsInitialValues();
}

void DriftModelSolver::Test()
{
	SIMPLE();
}

const std::valarray<double>& DriftModelSolver::GetV_m() const
{
	return _v_m;
}

const std::valarray<double>& DriftModelSolver::GetV_g() const
{
	return _v_g;
}

const std::valarray<double>& DriftModelSolver::GetV_l() const
{
	return _v_l;
}

const std::valarray<double>& DriftModelSolver::GetAlpha_g() const
{
	return _alpha_g;
}

const std::valarray<double>& DriftModelSolver::GetAlpha_l() const
{
	return _alpha_l;
}

const std::valarray<double>& DriftModelSolver::Get_P() const
{
	return _p;
}

const double & DriftModelSolver::GetDt() const
{
	return _dt;
}

const double & DriftModelSolver::GetDz() const
{
	return _dz;
}


void DriftModelSolver::InitializeGeometryParameters()
{

	size_t number_of_filled_cells = 0;
	size_t number_of_cells_in_current_segment = 0;
	size_t index_of_current_cell = 0;

	for (const WellSegment & wellSegment : _well.getSegments())
	{
		number_of_cells_in_current_segment = static_cast<int>(wellSegment.length / _dz);
		
		for (size_t i = 0; i < number_of_cells_in_current_segment; ++i)
		{

			index_of_current_cell = number_of_filled_cells + i;

			// ������������� ���������� �����
			_eps[index_of_current_cell] = wellSegment.relative_roughness;
			_d[index_of_current_cell] = wellSegment.diameter;
			_theta[index_of_current_cell] = wellSegment.tilt_angle;
		}

		number_of_filled_cells += number_of_cells_in_current_segment;
	}

}

void DriftModelSolver::initializeCellsInitialValues()
{
	std::cout << _n_points_cell_properties * _dz << std::endl;

	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{
		// ������ �������������� �������� ��� �������� ��������
		_alpha_g[i] = _drift_model.GetPRGasVolumeFractionInitialCondition();
		_alpha_l[i] = 1 - _alpha_g[i];
		_p[i] = i > 0 ? _p[i - 1] + _drift_model.CalculateHydrostaticPressure(_drift_model.GetMixtureDensity(_alpha_g[i], _alpha_l[i], 0), _dz) : 0.0;
	}

	

	for (size_t i = 0; i < _n_points_cell_velocities; ++i)
	{
		// ������ �������������� �������� ��� �������� ��������
		// �������� �������� ��� �������������� �������
		double v_steady = 0;
		_v_g[i] = 0;
		_v_l[i] = v_steady;
		_v_m[i] = _v_l[i];
	}
}


int DriftModelSolver::CalculateN()
{
	int n = 0;

	for (const WellSegment & wellSegment: _well.getSegments())
	{
		n += static_cast<int>(wellSegment.length / _dz);
	}

	return n;
}


std::valarray<double> DriftModelSolver::CalculateApproximateMixtureSpeed()
{
	// ��������
	const double accuracy = 0.01;
	// ������� ����������
	bool convergence_predicate;
	// ����������� ���������� (�������� ���. 106)
	double alpha_relax = 1;
	// �������� L2 ����� �������, ���������� �� �������� �������� �� ������� � ���������� ��������� �����
	double l2_norm_of_difference = 0.0;
	// �������� �������� ����� �� ���������� ��������� (� ������ ���������������� ��������� �� ���������� ��������� ����)
	std::valarray<double> v_m_star = alpha_relax * _v_m;
	// �������� �� ������� ��������
	std::valarray<double> v_m_current(_n_points_cell_velocities);
	// ������������ � �������
	std::valarray<double> alpha_p(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_e(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_w(0.0, _n_points_cell_velocities);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_velocities);


	do
	{
		for (size_t i = 0; i < _n_points_cell_velocities; i++)
		{
			// �������� � ����� P' (�������� � ��������� ���������)
			double p_P_stroke = _p[i];
			double alpha_g_P_stroke = _alpha_g[i];
			double alpha_l_P_stroke = _alpha_l[i];
			double rho_m_P_stroke = _drift_model.GetMixtureDensity(alpha_g_P_stroke, alpha_l_P_stroke, p_P_stroke);
			// �������� � ����� W' (�������� � ��������� ���������)
			double p_W_stroke = _p[i + 1];
			double alpha_g_W_stroke = _alpha_g[i + 1];
			double alpha_l_W_stroke = _alpha_l[i + 1];
			double rho_m_W_stroke = _drift_model.GetMixtureDensity(alpha_g_W_stroke, alpha_l_W_stroke, p_W_stroke);
			// �������� � ����� E' (�������� � ��������� ���������)
			double p_E_stroke = i > 0 ? _p[i - 1] : 0;
			double alpha_g_E_stroke = i > 0 ? _alpha_g[i - 1] : 0;
			double alpha_l_E_stroke = i > 0 ? _alpha_l[i - 1] : 0;
			double rho_m_E_stroke = _drift_model.GetMixtureDensity(alpha_g_E_stroke, alpha_l_E_stroke, p_E_stroke);
			// �������� � ����� P
			double v_m_star_P = v_m_star[i];
			double v_m_zero_P = _v_m[i];
			double p_P = (p_E_stroke + p_P_stroke) / 2;
			double alpha_g_P = _alpha_g[i];
			double alpha_l_P = _alpha_l[i];
			double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);
			double eps_P = _eps[i];
			double d_P = _d[i];
			double theta_P = _theta[i];
			double f_star_P = _drift_model.GetFrictionCoefficient(alpha_g_P, alpha_l_P, v_m_star_P, p_P, d_P, eps_P);
			// �������� � ����� W
			double v_m_star_W = i > 0 ? v_m_star[i - 1] : 0;
			// �������� � ����� E
			double v_m_star_E = i < _n_points_cell_velocities - 1 ? v_m_star[i + 1] : 0;
			// �������� �� ������ ����� (e) ��������� ������ ��� ����� P
			double v_m_star_e = (v_m_star_P + v_m_star_E) / 2;
			// �������� �� ����� ����� (w) ��������� ������ ��� ����� P
			double v_m_star_w = (v_m_star_W + v_m_star_P) / 2;

			alpha_e[i] = 0.5 * _dt / _dz * std::max(-v_m_star_e, 0.0);
			alpha_w[i] = 0.5 * _dt / _dz * std::max(v_m_star_w, 0.0);
			alpha_p[i] = 1 + 0.5 * _dt / _dz * (std::max(v_m_star_e, 0.0) +
						 std::max(-v_m_star_w, 0.0)) +
						 _dt * 2 * f_star_P * abs(v_m_star_P) / d_P;
			b[i] = v_m_zero_P
				+ _dt * _drift_model.g * cos(theta_P)
				- _dt / _dz * (2 / (rho_m_E_stroke + rho_m_P_stroke)) * (p_E_stroke - p_P_stroke)
				+ _dt / _dz * (2 / (rho_m_P_stroke + rho_m_W_stroke)) * (p_P_stroke - p_W_stroke);
		}
		

		TDMA(v_m_current, alpha_p, alpha_e, alpha_w, b);

		// �������� �������� �� ������� � ���������� ��������� �����
		l2_norm_of_difference = sqrt((v_m_star - v_m_current).apply([](double value)->double {return value * value; }).sum());
		convergence_predicate = l2_norm_of_difference > accuracy;
		v_m_star = v_m_current;

		// DEBUG
		// std::cout << "v_m approx" << std::endl << "\t\tv_m approx L2 norm of difference :" << l2_norm_of_difference << std::endl;

	} while (convergence_predicate);



	return v_m_current;
}

std::valarray<double> DriftModelSolver::CalculateMixtureVelocityCorrection(const std::valarray<double>& p_corr, const std::valarray<double> & v_m_approx)
{
	std::valarray<double> v_corr(0.0, _n_points_cell_velocities);

	for (size_t i = 0; i < _n_points_cell_velocities; ++i)
	{
		// �������� � ������� �����
		double p_P = p_corr[i];
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = _alpha_l[i];
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);

		// �������� � ����� ������ �� �������
		double p_E = p_corr[i + 1];
		double alpha_g_E = _alpha_g[i + 1];
		double alpha_l_E = _alpha_l[i + 1];
		double rho_m_E = _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E);

		
		// �������� ��������
		v_corr[i] = - 2 * ( _dt / (rho_m_P + rho_m_E)) * ((p_E - p_P) / _dz);
	}

	return v_corr;
}

std::valarray<double> DriftModelSolver::CalculateGasVelocity(const std::valarray<double> & p, const std::valarray<double> & v_m)
{
	std::valarray<double>v_g(0.0, _n_points_cell_velocities);
	double dp_dz;
	double v_d, C_0;
	double alpha_g_0 = 0;
	double p_wf = _drift_model.CalculateHydrostaticPressure(_drift_model.GetLiquidDensity(0), _well.GetLength());
	double U = _drift_model.GetPRCharacteristicVelocity(p_wf, _well.GetBottomCrossSectionArea());

	for (size_t i = 0; i < _n_points_cell_velocities; ++i)
	{
		dp_dz = (p[i + 1] - p[i]) / _dz;
		C_0 = _drift_model.CalculateGasProfileParameter(alpha_g_0, U, (_d[i] + _d[i + 1]) / 2, (p[i] + p[i + 1]) / 2, dp_dz);
		v_d = _drift_model.CalculateDriftVelocity(alpha_g_0, U, (_d[i] + _d[i + 1]) / 2, (p[i] + p[i + 1]) / 2);
		v_g[i] = C_0 * v_m[i] + v_d;
	}

	return v_g;
}

std::valarray<double> DriftModelSolver::CalculateLiquidVelocity(const std::valarray<double>& alpha_g, const std::valarray<double>& alpha_l, const std::valarray<double>& v_g, const std::valarray<double>& v_m)
{
	 std::valarray<double> v_l(_n_points_cell_velocities);

	 for (size_t i = 0; i < _n_points_cell_velocities; i++)
	 {
		 double alpha_gas = (alpha_g[i] + alpha_g[i + 1]) / 2;
		 double alpha_liquid = (alpha_l[i] + alpha_l[i + 1]) / 2;
		 v_l[i] = (v_m[i] - alpha_gas * v_g[i]) / alpha_liquid;
	 }

	 return v_l;
}

std::valarray<double> DriftModelSolver::CalculateGasVelocity_TEST(const std::valarray<double>& v_m)
{
	std::valarray<double>v_g(0.0, _n_points_cell_velocities);
	double v_d, C_0;

	for (size_t i = 0; i < _n_points_cell_velocities; ++i)
	{
		C_0 = _drift_model.CalculateGasProfileParameter_TEST(_alpha_g[i]);
		v_d = _drift_model.CalculateDriftVelocity_TEST(_alpha_g[i]);

		v_g[i] = C_0 * v_m[i] + v_d;
	}

	return v_g;
}


std::valarray<double> DriftModelSolver::CalculateGasVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_g)
{
	std::valarray<double>alpha_gas(0.0, _n_points_cell_properties);
	// ������������
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_properties);

	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{
		// �������� �� ������� ��������� ����     
		double alpha_g_0 = _alpha_g[i];
		double rho_g_0 = _drift_model.GetGasDensity(_p[i]);
		// �������� �� ������� ��������� ����   
		double rho_g = _drift_model.GetGasDensity(p_current[i]);

		// �������� �� ������ ����� ������������ ������     
		double v_g_e = i < _n_points_cell_properties - 1 ? v_g[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (p_current[i] + p_current[i + 1]) / 2 : 0;
		double rho_g_e = _drift_model.GetGasDensity(p_e);
		double F_e = rho_g_e * v_g_e;

		// �������� �� ����� ����� ���������� ������     
		double v_g_w = i > 0 ? v_g[i - 1] : 0;   
		double p_w = i > 0 ? (p_current[i] + p_current[i - 1]) / 2 : 0;
		double rho_g_w = _drift_model.GetGasDensity(p_w);
		double F_w = rho_g_w * v_g_w;							

		 alpha_e[i] = std::max(-F_e, 0.0);
		 alpha_w[i] = std::max(F_w, 0.0);
		 alpha_p[i] = rho_g * _dz / _dt + alpha_e[i] + alpha_w[i] + F_e - F_w;
		 b[i] = alpha_g_0 * rho_g_0 * _dz / _dt;
	}


	TDMA(alpha_gas, alpha_p, alpha_e, alpha_w, b);

	return alpha_gas;
}


std::valarray<double> DriftModelSolver::CalculateLiquidVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_l)
{
	std::valarray<double>alpha_liquid(0.0, _n_points_cell_properties);
	// ������������
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_properties);

	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{
		// �������� �� ������� ��������� ����     
		double alpha_l_0 = _alpha_l[i];
		double rho_l_0 = _drift_model.GetLiquidDensity(_p[i]);
		// �������� �� ������� ��������� ����   
		double rho_l = _drift_model.GetLiquidDensity(p_current[i]);

		// �������� �� ������ ����� ������������ ������     
		double v_l_e = i < _n_points_cell_properties - 1 ? v_l[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (p_current[i] + p_current[i + 1]) / 2 : 0;
		double rho_l_e = _drift_model.GetLiquidDensity(p_e);
		double F_e = rho_l_e * v_l_e;

		// �������� �� ����� ����� ���������� ������     
		double v_l_w = i > 0 ? v_l[i - 1] : 0;
		double p_w = i > 0 ? (p_current[i] + p_current[i - 1]) / 2 : 0;
		double rho_l_w = _drift_model.GetLiquidDensity(p_w);
		double F_w = rho_l_w * v_l_w;

		alpha_e[i] = std::max(-F_e, 0.0);
		alpha_w[i] = std::max(F_w, 0.0);
		alpha_p[i] = rho_l * _dz / _dt + alpha_e[i] + alpha_w[i] + F_e - F_w;
		b[i] = alpha_l_0 * rho_l_0 * _dz / _dt;
	}


	TDMA(alpha_liquid, alpha_p, alpha_e, alpha_w, b);

	return alpha_liquid;
}


std::valarray<double> DriftModelSolver::CalculatePressureCorrection(const std::valarray<double> & v_m_approx, const std::valarray<double> & alpha_g_past, const std::valarray<double> & alpha_l_past, const std::valarray<double> & p_past)
{
	// �������� �� ��������
	std::valarray<double> p_corr(0.0, _n_points_cell_properties);
	// ������������
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_properties);


	// ���������� ������������� 
	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{ 
		// �������� ��� ������� ��������� ������� � �������� ������
		double p_wf = _drift_model.CalculateHydrostaticPressure(_drift_model.GetLiquidDensity(0), _well.GetLength());
		double U = _drift_model.GetPRCharacteristicVelocity(p_wf, _well.GetBottomCrossSectionArea());

		/* ������� �������� ����� */

		// �������� � ������ ��������� ������ �� ���������� ��������� ����
		double alpha_g_0_P = alpha_g_past[i];
		double alpha_l_0_P = alpha_l_past[i];
		double p_past_0_P = p_past[i];
		double rho_m_0_P = _drift_model.GetMixtureDensity(alpha_g_0_P, alpha_l_0_P, p_past_0_P);

		// �������� �� ������� ��������� ����
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = _alpha_l[i];
		double p_P = _p[i];
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);

		// �������� �� ������ �����
		double v_m_star_e = i < _n_points_cell_velocities ? v_m_approx[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (_p[i] + _p[i + 1]) / 2 : 0;
		double alpha_g_e = i < _n_points_cell_properties - 1 ? (_alpha_g[i] + _alpha_g[i + 1]) / 2 : 0;
		double alpha_l_e = i < _n_points_cell_properties - 1 ? (_alpha_l[i] + _alpha_l[i + 1]) / 2 : 0;
		double rho_g_e = _drift_model.GetGasDensity(p_e);
		double rho_l_e = _drift_model.GetLiquidDensity(p_e);
		double rho_m_e = _drift_model.GetMixtureDensity(alpha_g_e, alpha_l_e, p_e);
		double d_e = i < _n_points_cell_properties - 1 ? (_d[i] + _d[i + 1]) / 2 : 1;
		double eps_e = i < _n_points_cell_properties - 1 ? (_eps[i] + _eps[i + 1]) / 2 : 0;
		double f_e = _drift_model.GetFrictionCoefficient(alpha_g_e, alpha_l_e, v_m_star_e, p_e, d_e, eps_e);
		
		// �������� ������� � �������� ������
		double dp_dz_e = i <_n_points_cell_properties - 1 ? (_p[i + 1] - _p[i]) / _dz : 0;
		double C_0_e = _drift_model.CalculateGasProfileParameter(alpha_g_e, U, d_e, p_e, dp_dz_e);
		double v_d_e = _drift_model.CalculateDriftVelocity(alpha_g_e, U, d_e, p_e);

		// �������� �� ����� ����� ������������ ������
		double v_m_star_w = i > 0 ? v_m_approx[i - 1] : 0;
		double p_w = i > 0 ? (_p[i] + _p[i - 1]) / 2 : 0;
		double alpha_g_w = i > 0 ? (_alpha_g[i] + _alpha_g[i - 1]) / 2 : 0;
		double alpha_l_w = i > 0 ? (_alpha_l[i] + _alpha_l[i - 1]) / 2 : 0;
		double rho_g_w = _drift_model.GetGasDensity(p_w);
		double rho_l_w = _drift_model.GetLiquidDensity(p_w);
		double rho_m_w = _drift_model.GetMixtureDensity(alpha_g_w, alpha_l_w, p_w);
		double d_w = i > 0 ? (_d[i] + _d[i - 1]) / 2 : 1;
		double eps_w = i > 0 ? (_eps[i] + _eps[i - 1]) / 2 : 0;
		double f_w = _drift_model.GetFrictionCoefficient(alpha_g_w, alpha_l_w, v_m_star_w, p_w, d_w, eps_w);

		// �������� ������� � �������� ������
		double dp_dz_w = i > 0  ? (_p[i] - _p[i - 1]) / _dz : 0;
		double C_0_w = _drift_model.CalculateGasProfileParameter(alpha_g_w, U, d_w, p_w, dp_dz_w);
		double v_d_w = _drift_model.CalculateDriftVelocity(alpha_g_w, U, d_w, p_w);


		/* �������� � ������ ��������� ������, ������������ ������ �� �������� */
		double alpha_g_E = i < _n_points_cell_properties - 1? _alpha_g[i + 1]: 0;
		double alpha_l_E = i < _n_points_cell_properties - 1 ? _alpha_l[i + 1] : 0;
		double p_E = i < _n_points_cell_properties - 1 ? _p[i + 1]: 0;
		double rho_m_E = i < _n_points_cell_properties - 1 ? _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E) : 0;
		


		/* �������� � ������ ��������� ������, ������������ ����� �� �������� */
		double alpha_g_W = i > 0 ? _alpha_g[i - 1] : 0;
		double alpha_l_W = i > 0 ? _alpha_l[i - 1] : 0;
		double p_W = i > 0 ? _p[i - 1] : 0;
		double rho_m_W = i > 0 ? _drift_model.GetMixtureDensity(alpha_g_W, alpha_l_W, p_W) : 0;
		

		// �������� �������� � ��������
		/*alpha_e[i] = (2 * _dt / (_dz * (rho_m_P + rho_m_E))) * (C_0_e * alpha_g_e * rho_g_e + (1 - alpha_g_e * C_0_e) * rho_l_e);
		alpha_w[i] = (2 * _dt / (_dz * (rho_m_W + rho_m_P))) * (C_0_w * alpha_g_w * rho_g_w + (1 - alpha_g_w * C_0_w) * rho_l_w);



		alpha_p[i] = alpha_e[i] + alpha_w[i];

		b[i] = _dz / _dt * (rho_m_0_P - rho_m_P)
			- v_d_e * alpha_g_e * rho_g_e
			- v_m_star_e * C_0_e * alpha_g_e * rho_g_e
			+ v_d_w * alpha_g_w * rho_g_w
			+ v_m_star_w * C_0_w * alpha_g_w * rho_g_w
			+ v_d_e * alpha_g_e * rho_l_e
			- v_m_star_e * ((1 - alpha_g_e * C_0_e)* rho_l_e)
			- v_d_w * alpha_g_w * rho_l_w
			+ v_m_star_w * ((1 - alpha_g_w * C_0_w)* rho_l_w);*/

		// ������ ���������� ���������
		alpha_e[i] = 2 * _dt / (_dz * (rho_m_P + rho_m_E));
		alpha_w[i] = 2 * _dt / (_dz * (rho_m_W + rho_m_P));
		alpha_p[i] = alpha_e[i] + alpha_w[i];
		b[i] = v_m_star_w - v_m_star_e;

	}

	TDMA(p_corr, alpha_p, alpha_e, alpha_w, b);

	return p_corr;
}

double DriftModelSolver::CalculateGasImbalance(const std::valarray<double>& alpha_g, const std::valarray<double>& alpha_g_past, const std::valarray<double>& p, const std::valarray<double>& p_past, const std::valarray<double>& v_g)
{
	double imbalance_value = 0;
	for (size_t i = 0; i < _n_points_cell_properties; i++)
	{
		// ������� �������� �����. cfv (current finite volume) 
		double v_w_cfv = i > 0 ? v_g[i - 1] : 0;
		double v_e_cfv = i < _n_points_cell_properties - 1 ? v_g[i] : 0;
		double v_cfv = (v_w_cfv + v_e_cfv) / 2;
		double alpha_cfv = alpha_g[i];
		double p_cfv = p[i];
		double rho_cfv = _drift_model.GetGasDensity(p_cfv);
		// ������� �������� ����� (�������� �� ���������� ��������� ����)
		double alpha_cfv_past = alpha_g_past[i];
		double p_cfv_past = p_past[i];
		double rho_cfv_past = _drift_model.GetGasDensity(p_cfv_past);
		// �������� ����� ����� �� ��������. lfv (left finite volume)
		double v_w_lfv = i > 1 ? v_g[i - 2] : 0;
		double v_e_lfv = (i > 0) && (i < _n_points_cell_properties - 1) ? v_g[i -1] : 0;
		double v_lfv = (v_w_lfv + v_e_lfv) / 2;
		double alpha_lfv = i > 0 ? alpha_g[i - 1] : 0;
		double p_lfv = i > 0 ? p[i - 1]: 0;
		double rho_lfv = _drift_model.GetGasDensity(p_lfv);

		// ��������� ��������� 
		imbalance_value += (alpha_cfv * rho_cfv - alpha_cfv_past * rho_cfv_past) / _dt + (alpha_cfv * rho_cfv * v_cfv - alpha_lfv * rho_lfv * v_lfv) / _dz;
	}

	return imbalance_value;
}

double DriftModelSolver::CalculateLiquidImbalance(const std::valarray<double>& alpha_l, const std::valarray<double>& alpha_l_past, const std::valarray<double>& p, const std::valarray<double>& p_past, const std::valarray<double>& v_l)
{
	double imbalance_value = 0;
	for (size_t i = 0; i < _n_points_cell_properties; i++)
	{
		// ������� �������� �����. cfv (current finite volume) 
		double v_w_cfv = i > 0 ? v_l[i - 1] : 0;
		double v_e_cfv = i < _n_points_cell_properties - 1 ? v_l[i] : 0;
		double v_cfv = (v_w_cfv + v_e_cfv) / 2;
		double alpha_cfv = alpha_l[i];
		double p_cfv = p[i];
		double rho_cfv = _drift_model.GetLiquidDensity(p_cfv);
		// ������� �������� ����� (�������� �� ���������� ��������� ����)
		double alpha_cfv_past = alpha_l_past[i];
		double p_cfv_past = p_past[i];
		double rho_cfv_past = _drift_model.GetLiquidDensity(p_cfv_past);
		// �������� ����� ����� �� ��������. lfv (left finite volume)
		double v_w_lfv = i > 1 ? v_l[i - 2] : 0;
		double v_e_lfv = (i > 0) && (i < _n_points_cell_properties - 1) ? v_l[i - 1] : 0;
		double v_lfv = (v_w_lfv + v_e_lfv) / 2;
		double alpha_lfv = i > 0 ? alpha_l[i - 1] : 0;
		double p_lfv = i > 0 ? p[i - 1] : 0;
		double rho_lfv = _drift_model.GetLiquidDensity(p_lfv);

		// ��������� ��������� 
		imbalance_value += (alpha_cfv * rho_cfv - alpha_cfv_past * rho_cfv_past) / _dt + (alpha_cfv * rho_cfv * v_cfv - alpha_lfv * rho_lfv * v_lfv) / _dz;
	}

	return imbalance_value;
}

// ����� ��������� ��������
void DriftModelSolver::TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d)
{
	// ������ ���������� �������
	size_t n = v.size();

	// ������������ � ������ ��������
	std::valarray<double> p(0.0, n), q(0.0, n);

	// ��������� �������� �������������
	p[0] = 0;
	q[0] = 0;
	

	// ���������� ������������� �� ����������� ��������
	for (size_t i = 1; i < n; ++i)
	{
		p[i] = b[i] / (a[i] - c[i] * p[i - 1]);
		q[i] = (c[i] * q[i - 1] + d[i]) / (a[i] - c[i] * p[i - 1]);
	}

	// �������� ���
	v[n - 1] = q[n - 1];

	for (int i = n - 2; i >= 0; --i)
	{
		v[i] = p[i] * v[i + 1] + q[i];
	}

}

void DriftModelSolver::SIMPLE()
{
	// ���������� ��� ������������� ����������
	std::valarray<double> v_m_intermediate(0.0, _n_points_cell_velocities);
	std::valarray<double> p_intermediate(0.0, _n_points_cell_properties);
	std::valarray<double> v_g_intermediate(0.0, _n_points_cell_velocities);
	std::valarray<double> v_l_intermediate(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_g_intermediate(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_l_intermediate(0.0, _n_points_cell_properties);

	// �������� � ������������� �������� �����
	std::valarray<double> v_m_approx;
	std::valarray<double> p_corr;
	std::valarray<double> v_m_corr;


	// �������� ����� �������� �������� ������� �� ������� � ���������� ��������� �����
	double l2_norm_of_difference_v_m = 0.0;
	double l2_norm_of_difference_p = 0.0;
	double l2_norm_of_difference_alpha_g = 0.0;

	// ��������
	const double accuracy = 0.1;
	bool norm_convergence_predicate;
	bool imbalance_convergence_predicate;


	// �������� (�������� 106) 
	double alpha_p_relax = 0.05;

	// ����� ���������� ��������
	int iteration_number = 0;

	// �������� � ���������� ��������
	std::valarray<double> p_past;
	std::valarray<double> alpha_g_past;
	std::valarray<double> alpha_l_past;

	// ����� �������� ��������
	double characteristic_stop_time = _drift_model.GetPRCharacteristicStopTime(_p[_n_points_cell_properties - 1], _well.GetBottomCrossSectionArea(), _well.GetLength());


	do
	{
		if ((iteration_number * _dt) >= characteristic_stop_time)
		{
			// �������� ��������
			_v_l[0] = 0.0;
			_v_g[0] = 0.0;
		}


		// ��������� �������
		_v_l[_n_points_cell_velocities - 1] = _drift_model.GetPRLiquidVelocityBoundaryCondition(_p[_n_points_cell_properties - 1], _well.GetBottomCrossSectionArea());
		_alpha_l[_n_points_cell_properties - 1] = _drift_model.GetPRLiquidVolumeFractionBoundaryCondition(_v_m[_n_points_cell_velocities - 1], _v_l[_n_points_cell_velocities - 1]);
		

		if (iteration_number == 0) {
			p_past = _p;
			alpha_g_past = _alpha_g;
			alpha_l_past = _alpha_l;
		}

		// ���������� ������������ �������� �������� �����
		v_m_approx =  CalculateApproximateMixtureSpeed();

		// ���������� �������� � ��������
		p_corr = CalculatePressureCorrection(v_m_approx, alpha_g_past, alpha_l_past, p_past);

		// � ��������� ������� ������������� �� ���������
		p_corr[0] = 0;
		p_corr[_n_points_cell_properties - 1] = 0;

		// ����������� ��������
		p_intermediate = _p + alpha_p_relax * p_corr;

		// ����������� ��������
		v_m_corr = CalculateMixtureVelocityCorrection(p_corr, v_m_approx);
		v_m_intermediate = v_m_approx + v_m_corr;

		// ���������� �������� ����
		v_g_intermediate = CalculateGasVelocity(p_intermediate, v_m_intermediate);

		// ���������� �������� ���� ����
		alpha_g_intermediate = CalculateGasVolumeFraction(p_intermediate, v_g_intermediate);

		// ���������� �������� ��������
		v_l_intermediate = CalculateLiquidVelocity(alpha_g_intermediate, _alpha_l, v_g_intermediate, v_m_intermediate);

		// ���������� �������� ���� ��������
		alpha_l_intermediate = CalculateLiquidVolumeFraction(p_intermediate, v_l_intermediate);

		// �������� �������� �� ������� � ���������� ��������� �����
		l2_norm_of_difference_v_m = sqrt((v_m_intermediate - _v_m).apply([](double value)->double {return value * value; }).sum());
		l2_norm_of_difference_p = sqrt((p_intermediate - _p).apply([](double value)->double {return value * value; }).sum());
		l2_norm_of_difference_alpha_g = sqrt((alpha_g_intermediate - _alpha_g).apply([](double value)->double {return value * value; }).sum());

		// �������� � ���������� ��������
		alpha_g_past = _alpha_g;
		alpha_l_past = _alpha_l;
		p_past = _p;

		// ���������� ���������� ������������� ����������
		_v_m = v_m_intermediate;
		_p = p_intermediate;
		_v_g = v_g_intermediate;
		_alpha_g = alpha_g_intermediate;
		_alpha_l = alpha_l_intermediate;

		// ���������
		double gas_imbalance_value = CalculateGasImbalance(_alpha_g, alpha_g_past, _p, p_past, _v_g);
		double liquid_imbalance_value = CalculateLiquidImbalance(_alpha_l, alpha_l_past, _p, p_past, _v_l);
		double imbalance_value = gas_imbalance_value + liquid_imbalance_value;

		//DEBUG
		std::cout << "v_m, p , alpha" << std::endl;
		std::cout << "\t\t iteration : " << iteration_number << std::endl;
		std::cout << "\t\t model time : " << _dt * iteration_number << " sec." << std::endl;
		std::cout << "\t\t imbalance value : " << imbalance_value << std::endl;
		std::cout << "\t\t v_m L2 norm of difference : " << l2_norm_of_difference_v_m << std::endl;
		std::cout << "\t\t p L2 norm of difference : " << l2_norm_of_difference_p << std::endl;
		std::cout << "\t\t alpha_g L2 norm of difference : " << l2_norm_of_difference_alpha_g << std::endl << std::endl;

		// �������� ����������
		norm_convergence_predicate = (l2_norm_of_difference_v_m > accuracy) || (l2_norm_of_difference_p > accuracy) || (l2_norm_of_difference_alpha_g > accuracy);
		imbalance_convergence_predicate = abs(imbalance_value) > accuracy;

		// ���������� ��������
		++iteration_number;
		//system("pause");

	} while (imbalance_convergence_predicate);
}



