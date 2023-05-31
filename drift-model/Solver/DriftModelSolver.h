#pragma once

#include <valarray>
#include <algorithm>

#include "..\Well\Well.h"
#include "..\MathModel\MathModel.h"


class DriftModelSolver
{
public:
	DriftModelSolver(double dz, double dt, const Well& well, MathModel::TaskType task_type);
	void Solve();
	const std::valarray<double> & GetV_m() const;
	const std::valarray<double> & GetV_g() const;
	const std::valarray<double> & GetV_l() const;
	const std::valarray<double> & GetAlpha_g() const;
	const std::valarray<double> & GetAlpha_l() const;
	const std::valarray<double> & Get_P() const;
	const double & GetDt() const;
	const double & GetDz() const;

private:

	const Well _well; // Скважина
	const double _dz; // Шаг по пространству
	const double _dt; // Шаг по времени

	MathModel::DriftModel _drift_model; // Класс математической модели
	size_t _n_points_cell_properties; // Число точек для значений всех параметров кроме скоростей
	size_t _n_points_cell_velocities; // Число точек для значений скоростей

	// Вектора с расчётными параметрами
	std::valarray<double> _theta; // Угол наклона трубы
	std::valarray<double> _d; // Диаметр трубы
	std::valarray<double> _eps; // Относительная шероховатость

	std::valarray<double> _v_m; // Скрость смеси
	std::valarray<double> _v_g; // Скрость газа
	std::valarray<double> _v_l; // Скрость жидкости
	std::valarray<double> _alpha_g; // Объёмная доля газа
	std::valarray<double> _alpha_l; // Объёмная доля жидкости
	std::valarray<double> _p; // Давление дисперсной среды


	void InitializeGeometryParameters(); // Инициализация параметров скважины
	
	int CalculateN();
	
	std::valarray<double> CalculateApproximateMixtureSpeed();
	
	std::valarray<double> CalculateMixtureVelocityCorrection(const std::valarray<double> & p_corr, const std::valarray<double> & v_m_approx);
	
	std::valarray<double> CalculateGasVelocity(const std::valarray<double> & p, const std::valarray<double> & v_m);
	std::valarray<double> CalculateLiquidVelocity(const std::valarray<double> & alpha_g, const std::valarray<double> & alpha_l, const std::valarray<double> & v_g, const std::valarray<double> & v_m);
	
	std::valarray<double> CalculateGasVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_g);
	std::valarray<double> CalculateLiquidVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_l);
	std::valarray<double> CalculatePressureCorrection(const std::valarray<double> & v_m_approx, const std::valarray<double> & alpha_g_past, const std::valarray<double> & alpha_l_past, const std::valarray<double> & p_past);

	// TEST
	std::valarray<double> CalculateGasVelocity_TEST(const std::valarray<double> & v_m);

	double CalculateGasImbalance(const std::valarray<double>& alpha_g, const std::valarray<double>& alpha_g_past, const std::valarray<double>& p_g, const std::valarray<double>& p_g_past, const std::valarray<double>& v_g);
	double CalculateLiquidImbalance(const std::valarray<double>& alpha_l, const std::valarray<double>& alpha_l_past, const std::valarray<double>& p_l, const std::valarray<double>& p_l_past, const std::valarray<double>& v_l);

	void TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d);
	void SIMPLE();
};

