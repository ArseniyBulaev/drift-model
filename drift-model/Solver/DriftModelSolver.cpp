#include "DriftModelSolver.h"

//DEBUG
#include <iostream>



DriftModelSolver::DriftModelSolver(double dz, double dt,  const Well & well, MathModel::TaskType task_type): _dz(dz), _dt(dt), _well(well), _drift_model(task_type, well)
{
	// Вычисление числа точек для значений параметров
	_n_points_cell_properties = CalculateNumberOfPoints();
	_n_points_cell_velocities = _n_points_cell_properties - 1;

	// Инициализация векторов
	
	// Параметры трубы
	_theta = std::valarray<double>(_n_points_cell_properties); // Угол наклона трубы
	_d = std::valarray<double>(_n_points_cell_properties); // Диаметр трубы
	_eps = std::valarray<double>(_n_points_cell_properties); // Относительная шероховатость
	// Расчётные параметры
	_v_m = std::valarray<double>(_n_points_cell_velocities); // Скрость смеси
	_v_g = std::valarray<double>(_n_points_cell_velocities); // Скрость газа
	_v_l = std::valarray<double>(_n_points_cell_velocities); // Скрость жидкости
	_alpha_g = std::valarray<double>(_n_points_cell_properties); // Объёмная доля газа
	_p = std::valarray<double>(_n_points_cell_properties); // Давление в дисперсной среде


	// Инициализация параметров ячеек
	InitializeGeometryParameters();
}

void DriftModelSolver::Solve()
{
	// Инициализация начальных условий
	_drift_model.SetInitialConditions(_alpha_g, _p, _v_m, _v_g, _v_l, _dz);

	// Решение задачи методом SIMPLE
	SimpleAlgorithm();
}


#pragma region To Delete
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



const std::valarray<double>& DriftModelSolver::Get_P() const
{
	return _p;
}

const double& DriftModelSolver::GetDt() const
{
	return _dt;
}

const double& DriftModelSolver::GetDz() const
{
	return _dz;
}
#pragma endregion

#pragma region Support

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

			// Инициализация параметров ячеек
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
void DriftModelSolver::CalculateApproximateMixtureSpeed(std::valarray<double>& v_m_intermediate)
{


	// Значения скорости смеси на предыдущей итерациии (Нужно для решения нелинейного уравнения)
	const std::valarray<double> & v_m_star = v_m_intermediate;

	// Коэффициенты в матрице
	std::valarray<double> alpha_p(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_e(0.0, _n_points_cell_velocities);
	std::valarray<double> alpha_w(0.0, _n_points_cell_velocities);
	// Правая часть
	std::valarray<double> b(0.0, _n_points_cell_velocities);



	for (size_t i = 0; i < _n_points_cell_velocities; i++)
	{
		// Значения в точке P' (давление и остальные параметры)
		double p_P_stroke = _p[i];
		double alpha_g_P_stroke = _alpha_g[i];
		double alpha_l_P_stroke = 1 - alpha_g_P_stroke;
		double rho_m_P_stroke = _drift_model.GetMixtureDensity(alpha_g_P_stroke, alpha_l_P_stroke, p_P_stroke);
		// Значения в точке W' (давление и остальные параметры)
		double p_W_stroke = _p[i + 1];
		double alpha_g_W_stroke = _alpha_g[i + 1];
		double alpha_l_W_stroke = 1 - alpha_g_W_stroke;
		double rho_m_W_stroke = _drift_model.GetMixtureDensity(alpha_g_W_stroke, alpha_l_W_stroke, p_W_stroke);
		// Значения в точке E' (давление и остальные параметры)
		double p_E_stroke = i > 0 ? _p[i - 1] : 0;
		double alpha_g_E_stroke = i > 0 ? _alpha_g[i - 1] : 0;
		double alpha_l_E_stroke = 1 - alpha_g_E_stroke;
		double rho_m_E_stroke = _drift_model.GetMixtureDensity(alpha_g_E_stroke, alpha_l_E_stroke, p_E_stroke);
		// Значения в точке P
		double v_m_star_P = v_m_star[i];
		double v_m_zero_P = _v_m[i];
		double p_P = (p_E_stroke + p_P_stroke) / 2;
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);
		double eps_P = _eps[i];
		double d_P = _d[i];
		double theta_P = _theta[i];
		double f_star_P = _drift_model.GetFrictionCoefficient(alpha_g_P, alpha_l_P, v_m_star_P, p_P, d_P, eps_P);
		// Значения в точке W
		double v_m_star_W = i > 0 ? v_m_star[i - 1] : 0;
		// Значения в точке E
		double v_m_star_E = i < _n_points_cell_velocities - 1 ? v_m_star[i + 1] : 0;
		// Значения на правой грани (e) конечного объёма для точки P
		double v_m_star_e = (v_m_star_P + v_m_star_E) / 2;
		// Значения на левой грани (w) конечного объёма для точки P
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

		// Релаксация для скорости (Фиолетовая книжка. страница 145)
		double alpha_u_relax = 0.6; // Коэффициент релаксации
		alpha_p[i] /= alpha_u_relax;
		b[i] += (1 - alpha_u_relax) * alpha_p[i] * v_m_intermediate[i];
	}

	// Нужно ли считать скорость в граничной ячейке ?
	alpha_w[_n_points_cell_velocities - 1] = 0;
	alpha_e[_n_points_cell_velocities - 1] = 0;
	alpha_p[_n_points_cell_velocities - 1] = 1;
	b[_n_points_cell_velocities - 1] = _v_m[_n_points_cell_velocities - 1];

	TDMA(v_m_intermediate, alpha_p, alpha_e, alpha_w, b);

}

std::valarray<double> DriftModelSolver::CalculateMixtureVelocityCorrection(const std::valarray<double>& p_corr)
{
	std::valarray<double> v_corr(0.0, _n_points_cell_velocities);

	for (size_t i = 0; i < _n_points_cell_velocities; ++i)
	{
		// Значения в текущей точке
		double p_P = p_corr[i];
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);

		// Значение в точке справа от текущей
		double p_E = p_corr[i + 1];
		double alpha_g_E = _alpha_g[i + 1];
		double alpha_l_E = 1 - alpha_g_E;
		double rho_m_E = _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E);


		// Линейная поправка
		v_corr[i] = -2 * (_dt / (rho_m_P + rho_m_E)) * ((p_E - p_P) / _dz);
	}

	// В граничной ячейке корректировка не требуется ?
	v_corr[_n_points_cell_velocities - 1] = 0;

	return v_corr;
}

std::valarray<double> DriftModelSolver::CalculateGasVolumeFraction(const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g)
{
	std::valarray<double>alpha_gas(0.0, _n_points_cell_properties);
	// Коэффициенты
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// Правая часть
	std::valarray<double> b(0.0, _n_points_cell_properties);

	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{
		// Значения на прошлом временном шаге     
		double alpha_g_0 = _alpha_g[i];
		double rho_g_0 = _drift_model.GetGasDensity(_p[i]);
		// Значения на текущем временном шаге   
		double rho_g = _drift_model.GetGasDensity(p_intermediate[i]);

		// Значения на правой грани контрольного объёма     
		double v_g_e = i < _n_points_cell_properties - 1 ? v_g[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (p_intermediate[i] + p_intermediate[i + 1]) / 2 : 0;
		double rho_g_e = _drift_model.GetGasDensity(p_e);
		double F_e = rho_g_e * v_g_e;

		// Значения на левой грани контрольно объёма     
		double v_g_w = i > 0 ? v_g[i - 1] : 0;
		double p_w = i > 0 ? (p_intermediate[i] + p_intermediate[i - 1]) / 2 : 0;
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

std::valarray<double> DriftModelSolver::CalculateGasVelocity(const std::valarray<double>& p, const std::valarray<double>& v_m)
{
	std::valarray<double>v_g;
	std::valarray<double> C_0 = _drift_model.CalculateC_0(_d, p, _n_points_cell_velocities, _dz);
	std::valarray<double> v_d = _drift_model.CalculateV_d(_d, p, _n_points_cell_velocities);
	v_g = C_0 * v_m + v_d;
	return v_g;
}

std::valarray<double> DriftModelSolver::CalculatePressureCorrection(const std::valarray<double>& v_m_intermediate)
{
	// Поправка на давление
	std::valarray<double> p_corr(0.0, _n_points_cell_properties);
	// Коэффициенты
	std::valarray<double> alpha_p(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_e(0.0, _n_points_cell_properties);
	std::valarray<double> alpha_w(0.0, _n_points_cell_properties);
	// Правая часть
	std::valarray<double> b(0.0, _n_points_cell_properties);

	// Параметр профиля и скорость дрейфа
	std::valarray<double> C_0 = _drift_model.CalculateC_0(_d, _p, _n_points_cell_velocities, _dz);
	std::valarray<double> v_d = _drift_model.CalculateV_d(_d, _p, _n_points_cell_velocities);

	// Заполнение коэффициентов 
	for (size_t i = 0; i < _n_points_cell_properties; ++i)
	{
		/* Текущей конечный объём */

		// Значения в центре конечного объёма на предыдущем временном шаге
		double alpha_g_0_P = _alpha_g[i];
		double alpha_l_0_P = 1 - alpha_g_0_P;
		double p_past_0_P = _p[i];
		double rho_m_0_P = _drift_model.GetMixtureDensity(alpha_g_0_P, alpha_l_0_P, p_past_0_P);

		// Значения на текущем временном шаге
		double alpha_g_P = _alpha_g[i];
		double alpha_l_P = 1 - alpha_g_P;
		double p_P = _p[i];
		double rho_m_P = _drift_model.GetMixtureDensity(alpha_g_P, alpha_l_P, p_P);

		// Значения на правой грани

		double v_m_star_e = i < _n_points_cell_velocities ? v_m_intermediate[i] : 0;
		double p_e = i < _n_points_cell_properties - 1 ? (_p[i] + _p[i + 1]) / 2 : 0;
		double alpha_g_e = i < _n_points_cell_properties - 1 ? (_alpha_g[i] + _alpha_g[i + 1]) / 2 : 0;
		double alpha_l_e = 1 - alpha_g_e;
		double rho_g_e = _drift_model.GetGasDensity(p_e);
		double rho_l_e = _drift_model.GetLiquidDensity(p_e);
		double rho_m_e = _drift_model.GetMixtureDensity(alpha_g_e, alpha_l_e, p_e);
		double d_e = i < _n_points_cell_properties - 1 ? (_d[i] + _d[i + 1]) / 2 : 1;
		double eps_e = i < _n_points_cell_properties - 1 ? (_eps[i] + _eps[i + 1]) / 2 : 0;
		double f_e = _drift_model.GetFrictionCoefficient(alpha_g_e, alpha_l_e, v_m_star_e, p_e, d_e, eps_e);

		// Параметр профиля и скорость дрейфа
		double C_0_e = i < _n_points_cell_velocities ? C_0[i] : 0;
		double v_d_e = i < _n_points_cell_velocities ? v_d[i] : 0;

		// Значения на левой грани контрольного объёма
		double v_m_star_w = i > 0 ? v_m_intermediate[i - 1] : 0;
		double p_w = i > 0 ? (_p[i] + _p[i - 1]) / 2 : 0;
		double alpha_g_w = i > 0 ? (_alpha_g[i] + _alpha_g[i - 1]) / 2 : 0;
		double alpha_l_w = 1 - alpha_g_w;
		double rho_g_w = _drift_model.GetGasDensity(p_w);
		double rho_l_w = _drift_model.GetLiquidDensity(p_w);
		double rho_m_w = _drift_model.GetMixtureDensity(alpha_g_w, alpha_l_w, p_w);
		double d_w = i > 0 ? (_d[i] + _d[i - 1]) / 2 : 1;
		double eps_w = i > 0 ? (_eps[i] + _eps[i - 1]) / 2 : 0;
		double f_w = _drift_model.GetFrictionCoefficient(alpha_g_w, alpha_l_w, v_m_star_w, p_w, d_w, eps_w);

		// Параметр профиля и скорость дрейфа
		double C_0_w = i > 0 ? C_0[i - 1] : 0;
		double v_d_w = i > 0 ? v_d[i - 1] : 0;


		/* Значения в центре конечного объёма, находящегося справа от текущего */
		double alpha_g_E = i < _n_points_cell_properties - 1 ? _alpha_g[i + 1] : 0;
		double alpha_l_E = 1 - alpha_g_E;
		double p_E = i < _n_points_cell_properties - 1 ? _p[i + 1] : 0;
		double rho_m_E = i < _n_points_cell_properties - 1 ? _drift_model.GetMixtureDensity(alpha_g_E, alpha_l_E, p_E) : 0;



		/* Значения в центре конечного объёма, находящегося слева от текущего */
		double alpha_g_W = i > 0 ? _alpha_g[i - 1] : 0;
		double alpha_l_W = 1 - alpha_g_W;
		double p_W = i > 0 ? _p[i - 1] : 0;
		double rho_m_W = i > 0 ? _drift_model.GetMixtureDensity(alpha_g_W, alpha_l_W, p_W) : 0;


		// Линейная поправка к давлению
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

		// Случай постоянной плотности
		alpha_e[i] = 2 * _dt / (_dz * (rho_m_P + rho_m_E));
		alpha_w[i] = 2 * _dt / (_dz * (rho_m_W + rho_m_P));
		alpha_p[i] = alpha_e[i] + alpha_w[i];
		b[i] = v_m_star_w - v_m_star_e;

	}

	TDMA(p_corr, alpha_p, alpha_e, alpha_w, b);

	// В граничной ячейке корректировка не требуется
	alpha_e[_n_points_cell_properties - 1] = 0;
	alpha_w[_n_points_cell_properties - 1] = 0;
	alpha_p[_n_points_cell_properties - 1] = 1;
	b[_n_points_cell_properties - 1] = 0;


	return p_corr;
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

double DriftModelSolver::CalculateGasImbalance(const std::valarray<double>& alpha_g_intermediate, const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g_intermediate)
{
	double imbalance_value = 0;
	for (size_t i = 0; i < _n_points_cell_properties; i++)
	{
		// Текущий конечный объём. cfv (current finite volume) 
		double v_w_cfv = i > 0 ? v_g_intermediate[i - 1] : 0;
		double v_e_cfv = i < _n_points_cell_properties - 1 ? v_g_intermediate[i] : 0;
		double v_cfv = (v_w_cfv + v_e_cfv) / 2;
		double alpha_cfv = alpha_g_intermediate[i];
		double p_cfv = p_intermediate[i];
		double rho_cfv = _drift_model.GetGasDensity(p_cfv);
		// Текущий конечный объём (Значения на предыдущем временном шаге)
		double alpha_cfv_past = _alpha_g[i];
		double p_cfv_past = _p[i];
		double rho_cfv_past = _drift_model.GetGasDensity(p_cfv_past);
		// Конечный объём слева от текущего. lfv (left finite volume)
		double v_w_lfv = i > 1 ? v_g_intermediate[i - 2] : 0;
		double v_e_lfv = (i > 0) && (i < _n_points_cell_properties - 1) ? v_g_intermediate[i - 1] : 0;
		double v_lfv = (v_w_lfv + v_e_lfv) / 2;
		double alpha_lfv = i > 0 ? alpha_g_intermediate[i - 1] : 0;
		double p_lfv = i > 0 ? p_intermediate[i - 1] : 0;
		double rho_lfv = _drift_model.GetGasDensity(p_lfv);

		// Суммарный дисбаланс 
		imbalance_value += (alpha_cfv * rho_cfv - alpha_cfv_past * rho_cfv_past) / _dt + (alpha_cfv * rho_cfv * v_cfv - alpha_lfv * rho_lfv * v_lfv) / _dz;
	}

	return imbalance_value;
}

// Метод матричной прогонки
void DriftModelSolver::TDMA(std::valarray<double>& v, const std::valarray<double>& a, const std::valarray<double>& b, const std::valarray<double>& c, const std::valarray<double>& d)
{
	// Размер расчётного вектора
	int n = static_cast<int>(v.size());

	// Коэффициенты в методе прогонки
	std::valarray<double> p(0.0, n), q(0.0, n);

	// Стартовые значения коэффициентов
	p[0] = 0;
	q[0] = 0;


	// Вычисление коэффициентов по рекурентным формулам
	for (size_t i = 1; i < n; ++i)
	{
		p[i] = b[i] / (a[i] - c[i] * p[i - 1]);
		q[i] = (c[i] * q[i - 1] + d[i]) / (a[i] - c[i] * p[i - 1]);
	}

	// Обратный ход
	v[n - 1] = q[n - 1];

	for (int i = n - 2; i >= 0; --i)
	{
		v[i] = p[i] * v[i + 1] + q[i];
	}

}


#pragma endregion



void DriftModelSolver::SimpleAlgorithm()
{

	// Поправки
	std::valarray<double> p_corr;
	std::valarray<double> v_m_corr;


	// Значения нормы разности векторов решений на текущем и предыдущем временных шагах
	double l2_norm_of_difference_of_alpha_g = 0.0;

	// Точность
	const double accuracy = 0.01;
	bool imbalance_convergence_predicate;


	// Патанкар (страница 106) 
	double alpha_p_relax = 0.1;

	// Номер внутренней итерации
	int internal_iteration_number = 0;
	int external_iteration_number = 0;

	// Промежуточные вычисления
	std::valarray<double> v_m_intermediate = _v_m;
	std::valarray<double> p_intermediate = _p;
	std::valarray<double> v_g_intermediate = _v_g;
	std::valarray<double> alpha_g_intermediate = _alpha_g;


	do
	{
		// Граничные условия
		_drift_model.SetBoundaryConditions(_alpha_g, _p, _v_m, _v_g, _v_l, _well, _dt);


		do
		{
			// Вычисление приближённого значения скорости смеси
			CalculateApproximateMixtureSpeed(v_m_intermediate);

			// Вычисление поправки к давлению
			p_corr = CalculatePressureCorrection(v_m_intermediate);

			// Исправление давления
			p_intermediate = p_intermediate + alpha_p_relax * p_corr;

			// Исправление скорости
			v_m_corr = CalculateMixtureVelocityCorrection(p_corr);
			v_m_intermediate = v_m_intermediate + v_m_corr;

			// Вычисление скорости газа
			v_g_intermediate = CalculateGasVelocity(p_intermediate, v_m_intermediate);

			// Вычисление объёмной доли газа
			alpha_g_intermediate = CalculateGasVolumeFraction(p_intermediate, v_g_intermediate);

			// Дисбаланс
			double imbalance_value = CalculateGasImbalance(alpha_g_intermediate, p_intermediate, v_g_intermediate);

			// INFO
			
			std::cout << "iteration : " << internal_iteration_number << std::endl;
			std::cout << "imbalance value : " << imbalance_value << std::endl;

			// Контроль сходимости
			imbalance_convergence_predicate = abs(imbalance_value) > accuracy;

			// Выполнение внутренней итерации
			++internal_iteration_number;

		} while (imbalance_convergence_predicate);
		
		internal_iteration_number = 0;
		// Выполнение внешней итерации (итерацци по времени)
		++external_iteration_number;

		// Вычисление нормы для проверки на сходимость
		l2_norm_of_difference_of_alpha_g = sqrt((_alpha_g - alpha_g_intermediate).apply([](double value)->double {return value * value; }).sum());
		

		_v_m = v_m_intermediate;
		_p = p_intermediate;
		_v_g = v_g_intermediate;
		_alpha_g = alpha_g_intermediate;
		
		

		// INFO
		std::cout << "\t\t alpha_g L2 norm of difference : " << l2_norm_of_difference_of_alpha_g << std::endl;
		std::cout << "\t\t model time : " << _dt * external_iteration_number << " sec." << std::endl << std::endl;

	} while (/*l2_norm_of_difference_of_alpha_g >= accuracy*/ _dt * external_iteration_number <= 1);

	
}



