#include "DriftModelSolver.h"


//DEBUG
#include <iostream>
#include <stdexcept>


DriftModelSolver::DriftModelSolver(double time,double dz, double dt,  const Well & well, MathModel::TaskType task_type): _calculation_time(time), _dz(dz), _dt(dt), _well(well), _drift_model(task_type, well)
{
	// ���������� ����� ����� ��� �������� ����������
	_n_points_cell_properties = CalculateNumberOfPoints();
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
	_p = std::valarray<double>(_n_points_cell_properties); // �������� � ���������� �����

	// ��������� ��� ������ � ����
	_print_step = _calculation_time ;
	_time_iterations_count = _calculation_time / _dt;
	_print_iteration_step = _print_step / _dt;
	
	// ������������� ���������� �����
	InitializeGeometryParameters();
}

void DriftModelSolver::Solve()
{
	// ������������� ��������� �������
	_drift_model.SetInitialConditions(_alpha_g, _p, _v_m, _v_g, _v_l, _theta, _dz);

	// ������� ������ ������� SIMPLE
	SimpleAlgorithm();
}

void DriftModelSolver::SimpleAlgorithm()
{
	// ��������
	std::valarray<double> p_corr;
	std::valarray<double> v_m_corr;

	// ��������
	const double accuracy = 1E-3;

	

	// �������� ����������
	bool imbalance_convergence_predicate;
	bool norm_convergence_predicate;

	// ����� ���������� ��������
	int internal_iteration_number = 0;
	int time_iteration_number = 0;

	_drift_model.SetBoundaryConditions(_alpha_g, _p, _v_m, _v_g, _v_l, _well, _dt);

	// ������������� ����������
	std::valarray<double> v_m_intermediate = _v_m;
	std::valarray<double> p_intermediate = _p;
	std::valarray<double> v_g_intermediate = _v_g;
	std::valarray<double> v_l_intermediate = _v_l;
	std::valarray<double> alpha_g_intermediate = _alpha_g;

	// �������� � ���������� �������� ��� ������ ����������
	std::valarray<double> v_m_previous_iteration = _v_m;
	std::valarray<double> p_previous_iteration = _p;
	std::valarray<double> v_g_previous_iteration = _v_g;
	std::valarray<double> v_l_previous_iteration = _v_l;
	std::valarray<double> alpha_g_previous_iteration = _alpha_g;

	double l2_norm_of_difference_v_m;
	double l2_norm_of_difference_p;
	double l2_norm_of_difference_alpha_g;

	double l2_norm_of_v_m_previous_iteration;
	double l2_norm_of_p_previous_iteration;
	double l2_norm_of_alpha_g_previous_iteration;

	/*_results_writer.WriteToFile((_p) / _drift_model.atm, _dz, "p__" + std::to_string(time_iteration_number * _dt) + ".txt");
	_results_writer.WriteToFile(_alpha_g, _dz, "alpha_g__" + std::to_string(time_iteration_number * _dt) + ".txt");
	_results_writer.WriteToFile(_v_m, _dz, "v_m__" + std::to_string(time_iteration_number * _dt) + ".txt");*/


	do
	{
		std::cout << "\t\t model time : " << _dt * time_iteration_number << " sec." << std::endl << std::endl;

		// ��������� �������
		_drift_model.SetBoundaryConditions(_alpha_g, _p, _v_m, _v_g, _v_l, _well, _dt);

		do
		{
			// ������������� ���� �� ������� �������� ������� �������
 			// PrintCourantNumber(v_m_intermediate);

			// ���������� ������������ �������� �������� �����
			CalculateApproximateMixtureVelocity(v_m_intermediate);

			// ���������� �������� � ��������
			p_corr = CalculatePressureCorrection(v_m_intermediate, p_intermediate);

			// ����������� ��������
			p_intermediate = _p + alpha_p_relax * p_corr;

			// ����������� ��������
			v_m_corr = CalculateMixtureVelocityCorrection(p_corr);
			v_m_intermediate = v_m_intermediate + v_m_corr;

			// ���������� �������� ����
			v_g_intermediate = CalculateGasVelocity(p_intermediate, v_m_intermediate, alpha_g_previous_iteration);

			// ���������� �������� ���� ����
			alpha_g_intermediate = CalculateGasVolumeFraction(p_intermediate, v_g_intermediate);

			// ���������� �������� ��������
			v_l_intermediate = CalculateLiquidVelocity(v_m_intermediate, alpha_g_intermediate, v_g_intermediate);

			// ���������
			double imbalance_value = CalculateGasImbalance(alpha_g_intermediate, p_intermediate, v_g_intermediate, v_l_intermediate);
			imbalance_convergence_predicate = imbalance_value > accuracy;
			std::cout << "imbalance value : " << imbalance_value << std::endl;
			std::cout << "iteration : " << internal_iteration_number << std::endl;

			// ���������� ����� �������� �������
			l2_norm_of_difference_v_m = sqrt((v_m_intermediate - v_m_previous_iteration).apply([](double value)->double {return value * value; }).sum());
			l2_norm_of_difference_p = sqrt((p_intermediate - p_previous_iteration).apply([](double value)->double {return value * value; }).sum());
			l2_norm_of_difference_alpha_g = sqrt((alpha_g_intermediate - alpha_g_previous_iteration).apply([](double value)->double {return value * value; }).sum());
			// double l2_norm_of_difference_v_g = sqrt((v_g_intermediate - v_g_previous_iteration).apply([](double value)->double {return value * value; }).sum());

			// ���������� ����� ������� � ���������� ��������
			l2_norm_of_v_m_previous_iteration = sqrt((v_m_previous_iteration).apply([](double value)->double {return value * value; }).sum());
			l2_norm_of_p_previous_iteration = sqrt((p_previous_iteration).apply([](double value)->double {return value * value; }).sum());
			l2_norm_of_alpha_g_previous_iteration = sqrt((alpha_g_previous_iteration).apply([](double value)->double {return value * value; }).sum());

			// INFO
			std::cout << "\t\t v_m L2 norm of difference : " << l2_norm_of_difference_v_m / l2_norm_of_v_m_previous_iteration << std::endl;
			std::cout << "\t\t p L2 norm of difference : " << l2_norm_of_difference_p / l2_norm_of_p_previous_iteration << std::endl;
			//std::cout << "\t\t v_g L2 norm of difference : " << l2_norm_of_difference_v_g << std::endl;
			std::cout << "\t\t alpha_g L2 norm of difference : " << l2_norm_of_difference_alpha_g / l2_norm_of_alpha_g_previous_iteration << std::endl << std::endl;
			
			
			// �������� ���������� �� �����
			norm_convergence_predicate = ((l2_norm_of_difference_v_m / l2_norm_of_v_m_previous_iteration) > accuracy) ||
				 ((l2_norm_of_difference_p / l2_norm_of_p_previous_iteration) >  accuracy) ||
				 ((l2_norm_of_difference_alpha_g / l2_norm_of_alpha_g_previous_iteration) > accuracy);

			// ���������� ���������� ��������
			++internal_iteration_number;
		
			// ���������� �������� � ���������� ��������
			v_m_previous_iteration = v_m_intermediate;
			v_g_previous_iteration = v_g_intermediate;
			v_l_previous_iteration = v_l_intermediate;
			p_previous_iteration = p_intermediate;
			alpha_g_previous_iteration = alpha_g_intermediate;

		} while (norm_convergence_predicate);

		internal_iteration_number = 0;

		_v_m = v_m_intermediate;
		_p = p_intermediate;
		_v_g = v_g_intermediate;
		_v_l = v_l_intermediate;
		_alpha_g = alpha_g_intermediate;

		// ��������� �� ��������� ������ ������
		if (abs(_alpha_g[0]) > 0)
		{
			_alpha_g[0] += _alpha_g[0] * _v_g[0] * _dt / _dz;
		}

		// ���������� ������� �������� (�������� �� �������)
		++time_iteration_number;

		if (time_iteration_number % _print_iteration_step == 0)
		{
			

			_results_writer.WriteToFile((_p) / _drift_model.atm, _dz, "p__" + std::to_string(_dz) + ".txt");
			_results_writer.WriteToFile(_alpha_g, _dz, "alpha_g__" + std::to_string(_dz) + ".txt");
			_results_writer.WriteToFile(_v_m, _dz, "v_m__" + std::to_string(_dz) + ".txt");
			//_results_writer.WriteToFile(_v_g, _dz, "v_g__" + std::to_string(_dz) + ".txt");
			//_results_writer.WriteToFile(_v_l, _dz, "v_l__" + std::to_string(_dz) + ".txt");
		}
		

	} while (_dt * time_iteration_number <= _calculation_time);
}

#pragma region Support

void DriftModelSolver::PrintCourantNumber(const std::valarray<double>  & v_m_intermediate)
{
	// ���������� ����� �������
	auto U = v_m_intermediate.max();
	double C = abs(U) * _dt / _dz;
	std::cout << "Courant number value : " << C << std::endl;
	if (C > 1)
	{
		system("pause");
	}

}

void DriftModelSolver::InitializeGeometryParameters()
{

	size_t number_of_filled_cells = 0;
	size_t number_of_cells_in_current_segment = 0;
	size_t index_of_current_cell = 0;

	for (const WellSegment& wellSegment : _well.getSegments())
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

int DriftModelSolver::CalculateNumberOfPoints()
{
	int n = 0;

	for (const WellSegment& wellSegment : _well.getSegments())
	{
		n += static_cast<int>(wellSegment.length / _dz);
	}

	return n;
}
#pragma endregion

#pragma region Calculation

#pragma region Numeric method
void DriftModelSolver::CalculateApproximateMixtureVelocity(std::valarray<double>& v_m_intermediate)
{
	// �������� �������� ����� �� ���������� ��������� (����� ��� ������� ����������� ���������)
	const std::valarray<double>& v_m_star = v_m_intermediate;

	// ������������ � �������
	std::valarray<double> alpha_p(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_e(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_w(0.0, _n_points_cell_velocities);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_velocities);


	for (int i = 0; i < _n_points_cell_velocities; i++)
	{
		// �������� � ����� P' (�������� � ��������� ���������)
		double p_P_stroke = _p[i]; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_g_P_stroke = _alpha_g[i]; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_l_P_stroke = 1 - alpha_g_P_stroke;
		double rho_m_P_stroke = _drift_model.GetMixtureDensity(alpha_g_P_stroke, alpha_l_P_stroke, p_P_stroke, i * _dz);
		// �������� � ����� E' (�������� � ��������� ���������)
		double p_E_stroke = _p[i + 1]; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_g_E_stroke = _alpha_g[i + 1]; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_l_E_stroke = 1 - alpha_g_E_stroke;
		double rho_m_E_stroke = _drift_model.GetMixtureDensity(alpha_g_E_stroke, alpha_l_E_stroke, p_E_stroke, i * _dz);
		// �������� � ����� W' (�������� � ��������� ���������)
		double p_W_stroke = i > 0 ? _p[i - 1] : 0; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_g_W_stroke = i > 0 ? _alpha_g[i - 1] : 0; // ����� � ����������� ���������� ���� ��� � ���������� �������� ?
		double alpha_l_W_stroke = 1 - alpha_g_W_stroke;
		double rho_m_W_stroke = _drift_model.GetMixtureDensity(alpha_g_W_stroke, alpha_l_W_stroke, p_W_stroke, i * _dz);
		// �������� � ����� P
		double v_m_star_P = v_m_star[i];
		double v_m_zero_P = _v_m[i];
		double p_P = (p_E_stroke + p_P_stroke) / 2;
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P, i * _dz);
		double eps_P = (_eps[i] + _eps[i + 1]) / 2;
		double d_P = (_d[i] + _d[i + 1]) / 2;
		double theta_P = (_theta[i] + _theta[i + 1]) / 2;
		double f_star_P = _drift_model.GetFrictionCoefficient(alpha_g_P, alpha_l_P, v_m_star_P, p_P, d_P, eps_P);
		// �������� � ����� W
		double v_m_star_W = i > 0 ? v_m_star[i - 1] : v_m_star[i];
		// �������� � ����� E
		double v_m_star_E = i < _n_points_cell_velocities - 1 ? v_m_star[i + 1] : v_m_star[i];
		// �������� �� ������ ����� (e) ��������� ������ ��� ����� P
		double v_m_star_e = (v_m_star_P + v_m_star_E) / 2;
		// �������� �� ����� ����� (w) ��������� ������ ��� ����� P
		double v_m_star_w = (v_m_star_W + v_m_star_P) / 2;

		alpha_e[i] = 0.5 * _dt / _dz * std::max(-rho_m_E_stroke * v_m_star_e, 0.0);
		alpha_w[i] = 0.5 * _dt / _dz * std::max(rho_m_W_stroke*v_m_star_w, 0.0);
		alpha_p[i] = rho_m_P_stroke + 0.5 * _dt / _dz * (std::max(rho_m_E_stroke * v_m_star_e, 0.0) +
			std::max(-rho_m_W_stroke*v_m_star_w, 0.0)) +
			_dt * 2 * f_star_P * rho_m_P_stroke * abs(v_m_star_P) / d_P;


		b[i] = v_m_zero_P * rho_m_P_stroke
			+ (_dt * rho_m_P_stroke * _drift_model.g * cos(theta_P))
			- _dt / _dz * (p_E_stroke - p_P_stroke);

		// ���������� ��� �������� (���������� ������. �������� 145)
		alpha_p[i] /= alpha_v_relax;
		b[i] += (1 - alpha_v_relax) * alpha_p[i] * _v_m[i];
	}

	alpha_e[_n_points_cell_velocities - 1] = 0;
	alpha_w[_n_points_cell_velocities - 1] = 0;
	alpha_p[_n_points_cell_velocities - 1] = 1;
	b[_n_points_cell_velocities - 1] = _v_m[_n_points_cell_velocities - 1];

	TDMA(v_m_intermediate, alpha_p, alpha_e, alpha_w, b);
}
std::valarray<double> DriftModelSolver::CalculateMixtureVelocityCorrection(const std::valarray<double>& p_corr)
{
	std::valarray<double> v_corr(0.0, _n_points_cell_velocities);

	for (int i = 0; i < _n_points_cell_velocities; ++i)
	{
		// �������� � ������� �����
		double p_P = p_corr[i];
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P, _dz * i);

		// �������� � ����� ������ �� �������
		double p_E = p_corr[i + 1];
		double alpha_g_E = _alpha_g[i + 1];
		double alpha_l_E = 1 - alpha_g_E;
		double rho_m_E = _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E, _dz * i);


		// �������� ��������
		v_corr[i] = -2 * (_dt / (rho_m_P + rho_m_E)) * ((p_E - p_P) / _dz);
	}

	v_corr[_n_points_cell_velocities - 1] = 0;

	return v_corr;
}
std::valarray<double> DriftModelSolver::CalculateGasVolumeFraction(const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g)
{
	std::valarray<double>alpha_gas(0.0, _n_points_cell_properties);
	// ������������
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_properties);

	for (int i = 0; i < _n_points_cell_properties; ++i)
	{
		// �������� �� ������� ��������� ����     
		double alpha_g_0 = _alpha_g[i];
		double rho_g_0 = _drift_model.GetGasDensity(_p[i], _dz * i);
		// �������� �� ������� ��������� ����   
		double rho_g = _drift_model.GetGasDensity(p_intermediate[i], _dz * i);

		// �������� �� ������ ����� ������������ ������     
		double v_g_e = i < _n_points_cell_properties - 1 ? v_g[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (p_intermediate[i] + p_intermediate[i + 1]) / 2 : 0;
		double rho_g_e = _drift_model.GetGasDensity(p_e, _dz * i);
		double F_e = rho_g_e * v_g_e;

		// �������� �� ����� ����� ���������� ������     
		double v_g_w = i > 0 ? v_g[i - 1] : 0;
		double p_w = i > 0 ? (p_intermediate[i] + p_intermediate[i - 1]) / 2 : 0;
		double rho_g_w = _drift_model.GetGasDensity(p_w, _dz * i);
		double F_w = rho_g_w * v_g_w;

		alpha_e[i] = std::max(-F_e, 0.0);
		alpha_w[i] = std::max(F_w, 0.0);
		alpha_p[i] = rho_g * _dz / _dt + alpha_e[i] + alpha_w[i] + F_e - F_w;
		b[i] = alpha_g_0 * rho_g_0 * _dz / _dt;
	}


	alpha_e[_n_points_cell_properties - 1] = 0;
	alpha_w[_n_points_cell_properties - 1] = 0;
	alpha_p[_n_points_cell_properties - 1] = 1;
	b[_n_points_cell_properties - 1] = _alpha_g[_n_points_cell_properties - 1];

	TDMA(alpha_gas, alpha_p, alpha_e, alpha_w, b);

	// �������� �� ����� �� ������� ���������� ��������
	alpha_gas.apply([](double value)->double {

		if (value > 1)
		{
			value = 1;
		}

		if (value < 0)
		{
			value = 0;
		}

		return value; 
	});

	return alpha_gas;
}
std::valarray<double> DriftModelSolver::CalculatePressureCorrection(const std::valarray<double>& v_m_intermediate, const std::valarray<double> & p_intermediate)
{
	// �������� �� ��������
	std::valarray<double> p_corr(0.0, _n_points_cell_properties);
	// ������������
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// ������ �����
	std::valarray<double> b(0.0, _n_points_cell_properties);

	// �������� ������� � �������� ������
	std::valarray<double> C_0 = _drift_model.CalculateC_0(_d, _alpha_g, _p, _n_points_cell_velocities, _dz);
	std::valarray<double> v_d = _drift_model.CalculateV_d(_d, _alpha_g, _p, _n_points_cell_velocities);

	// ���������� ������������� 
	for (int i = 0; i < _n_points_cell_properties; ++i)
	{
		/* ������� �������� ����� */

		// �������� � ������ ��������� ������ �� ���������� ��������� ����
		double alpha_g_0_P = _alpha_g[i];
		double alpha_l_0_P = 1 - alpha_g_0_P;
		double p_past_0_P = _p[i];
		double rho_m_0_P = _drift_model.GetMixtureDensity(alpha_g_0_P, alpha_l_0_P, p_past_0_P, _dz * i);

		// �������� �� ������� ��������� ����
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double p_P = p_intermediate[i];
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P, _dz * i);

		// �������� �� ������ �����

		double v_m_star_e = i < _n_points_cell_velocities ? v_m_intermediate[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (p_intermediate[i] + p_intermediate[i + 1]) / 2 : 0;
		double alpha_g_e = i < _n_points_cell_properties - 1 ? (_alpha_g[i] + _alpha_g[i + 1]) / 2 : 0;
		double alpha_l_e = 1 - alpha_g_e;
		double rho_g_e = _drift_model.GetGasDensity(p_e, _dz * i);
		double rho_l_e = _drift_model.GetLiquidDensity(p_e);
		double rho_m_e = _drift_model.GetMixtureDensity(alpha_g_e, alpha_l_e, p_e, _dz * i);
		double d_e = i < _n_points_cell_properties - 1 ? (_d[i] + _d[i + 1]) / 2 : 1;
		double eps_e = i < _n_points_cell_properties - 1 ? (_eps[i] + _eps[i + 1]) / 2 : 0;
		double f_e = _drift_model.GetFrictionCoefficient(alpha_g_e, alpha_l_e, v_m_star_e, p_e, d_e, eps_e);

		// �������� ������� � �������� ������
		double C_0_e = i < _n_points_cell_velocities ? C_0[i] : 0;
		double v_d_e = i < _n_points_cell_velocities ? v_d[i] : 0;

		// �������� �� ����� ����� ������������ ������
		double v_m_star_w = i > 0 ? v_m_intermediate[i - 1] : 0;
		double p_w = i > 0 ? (p_intermediate[i] + p_intermediate[i - 1]) / 2 : 0;
		double alpha_g_w = i > 0 ? (_alpha_g[i] + _alpha_g[i - 1]) / 2 : 0;
		double alpha_l_w = 1 - alpha_g_w;
		double rho_g_w = _drift_model.GetGasDensity(p_w, _dz * i);
		double rho_l_w = _drift_model.GetLiquidDensity(p_w);
		double rho_m_w = _drift_model.GetMixtureDensity(alpha_g_w, alpha_l_w, p_w, _dz * i);
		double d_w = i > 0 ? (_d[i] + _d[i - 1]) / 2 : 1;
		double eps_w = i > 0 ? (_eps[i] + _eps[i - 1]) / 2 : 0;
		double f_w = _drift_model.GetFrictionCoefficient(alpha_g_w, alpha_l_w, v_m_star_w, p_w, d_w, eps_w);

		// �������� ������� � �������� ������
		double C_0_w = i > 0 ? C_0[i - 1] : 0;
		double v_d_w = i > 0 ? v_d[i - 1] : 0;


		/* �������� � ������ ��������� ������, ������������ ������ �� �������� */
		double alpha_g_E = i < _n_points_cell_properties - 1 ? _alpha_g[i + 1] : 0;
		double alpha_l_E = 1 - alpha_g_E;
		double p_E = i < _n_points_cell_properties - 1 ? p_intermediate[i + 1] : 0;
		double rho_m_E = i < _n_points_cell_properties - 1 ? _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E, _dz * i) : 0;



		/* �������� � ������ ��������� ������, ������������ ����� �� �������� */
		double alpha_g_W = i > 0 ? _alpha_g[i - 1] : 0;
		double alpha_l_W = 1 - alpha_g_W;
		double p_W = i > 0 ? p_intermediate[i - 1] : 0;
		double rho_m_W = i > 0 ? _drift_model.GetMixtureDensity(alpha_g_W, alpha_l_W, p_W, _dz * i) : 0;


		// �������� �������� � ��������
		alpha_e[i] = (2 * _dt / (_dz * (rho_m_P + rho_m_E))) * (C_0_e * alpha_g_e * rho_g_e + (1 - alpha_g_e * C_0_e) * rho_l_e);
		alpha_w[i] = (2 * _dt / (_dz * (rho_m_W + rho_m_P))) * (C_0_w * alpha_g_w * rho_g_w + (1 - alpha_g_w * C_0_w) * rho_l_w);
		alpha_p[i] = alpha_e[i] + alpha_w[i];

		b[i] = _dz / _dt * (rho_m_0_P - rho_m_P)
			- v_d_e * alpha_g_e * rho_g_e
			- v_m_star_e * C_0_e * alpha_g_e * rho_g_e
			+ v_d_w * alpha_g_w * rho_g_w
			+ v_m_star_w * C_0_w * alpha_g_w * rho_g_w
			+ v_d_e * alpha_g_e * rho_l_e
			- v_m_star_e * ((1 - alpha_g_e * C_0_e) * rho_l_e)
			- v_d_w * alpha_g_w * rho_l_w
			+ v_m_star_w * ((1 - alpha_g_w * C_0_w) * rho_l_w);
	}

	
	alpha_e[0] = 0;
	alpha_w[0] = 0;
	alpha_p[0] = 1;
	b[0] = 0; 

	alpha_e[_n_points_cell_properties - 1] = 0;
	alpha_w[_n_points_cell_properties - 1] = 0;
	alpha_p[_n_points_cell_properties - 1] = 1;
	b[_n_points_cell_properties - 1] = 0;

	

	TDMA(p_corr, alpha_p, alpha_e, alpha_w, b);

	return p_corr;
}
#pragma endregion

#pragma region Explicit
std::valarray<double> DriftModelSolver::CalculateGasVelocity(const std::valarray<double>& p, const std::valarray<double>& v_m, const std::valarray<double> & alpha_g_previous_iteration)
{
	std::valarray<double>v_g;
	std::valarray<double> C_0 = _drift_model.CalculateC_0(_d, alpha_g_previous_iteration,  p, _n_points_cell_velocities, _dz);
	std::valarray<double> v_d = _drift_model.CalculateV_d(_d, alpha_g_previous_iteration,  p, _n_points_cell_velocities);
	v_g = C_0 *v_m + v_d;
	return v_g;
}
std::valarray<double> DriftModelSolver::CalculateLiquidVelocity(const std::valarray<double>&v_m, const std::valarray<double>&alpha_g, const std::valarray<double>&v_g)
{
	std::valarray<double>v_l(_n_points_cell_velocities);

	for (size_t i = 0; i < _n_points_cell_velocities; i++)
	{
		double alpha_g_mid = (alpha_g[i] + alpha_g[i + 1]) / 2;
		double alpha_l_mid = 1 - alpha_g_mid;

		if (alpha_l_mid == 0)
		{
			throw std::invalid_argument("alpha_l can't be zero");
		}

		v_l[i] = (v_m[i] - alpha_g_mid * v_g[i]) / alpha_l_mid;
	}

	return v_l;
}
#pragma endregion

#pragma region Support
double DriftModelSolver::CalculateGasImbalance(const std::valarray<double>& alpha_g_intermediate, const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g_intermediate, const std::valarray<double>& v_l_intermediate)
{
	double gas_imbalance_value = 0;

	for (int i = 1; i < _n_points_cell_properties; i++)
	{
		// ����� P
		double rho_g_P = _drift_model.GetGasDensity(p_intermediate[i], _dz * i);
		double alpha_g_P = alpha_g_intermediate[i];
		double v_g_P = (v_g_intermediate[i - 1] + (i < _n_points_cell_properties - 1 ? v_g_intermediate[i] : 0)) / 2;
		double rho_g_o_P = _drift_model.GetGasDensity(_p[i], _dz * i);
		double alpha_g_o_P = _alpha_g[i];

		// ����� W
		double rho_g_W = _drift_model.GetGasDensity(p_intermediate[i - 1], _dz * i);
		double alpha_g_W = alpha_g_intermediate[i - 1];
		double v_g_W = ((i - 2 > 0 ? v_g_intermediate[i - 2] : 0) + v_g_intermediate[i - 1]) / 2;
	}

	return gas_imbalance_value;
}
// ����� ��������� ��������
void DriftModelSolver::TDMA(std::valarray<double>& v, const std::valarray<double>& a, const std::valarray<double>& b, const std::valarray<double>& c, const std::valarray<double>& d)
{
	// ������ ���������� �������
	int n = static_cast<int>(v.size());

	// ������������ � ������ ��������
	std::valarray<double> p(0.0, n), q(0.0, n);

	// ��������� �������� �������������
	p[0] = b[0] / a[0];
	q[0] = d[0] / a[0];


	// ���������� ������������� �� ����������� ��������
	for (int i = 1; i < n; ++i)
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
#pragma endregion

#pragma endregion
