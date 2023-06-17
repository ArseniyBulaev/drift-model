#pragma once

#include <valarray>
#include <algorithm>

#include "..\Well\Well.h"
#include "..\MathModel\MathModel.h"
#include"..\Writer.h"


class DriftModelSolver
{
public:
	DriftModelSolver(double time, double dz, double dt, const Well& well, MathModel::TaskType task_type);
	void Solve();

private:

	const Well _well; // Скважина
	const double _dz; // Шаг по пространству
	double _dt; // Шаг по времени
	double _calculation_time; // Время расчёта
	double _print_step; // Шаг печати в файл
	int _time_iterations_count;
	int _print_iteration_step;

	MathModel::DriftModel _drift_model; // Класс математической модели
	Writer _results_writer; // Объект для печати результатов

	size_t _n_points_cell_properties; // Число точек для значений всех параметров кроме скоростей
	size_t _n_points_cell_velocities; // Число точек для значений скоростей

	double alpha_p_relax = 0.9; // Коэффициент релаксации для давления
	double alpha_v_relax = 0.9; // Коэффициент релаксации для скорости

	// Вектора с расчётными параметрами
	std::valarray<double> _theta; // Угол наклона трубы
	std::valarray<double> _d; // Диаметр трубы
	std::valarray<double> _eps; // Относительная шероховатость

	std::valarray<double> _v_m; // Скрость смеси
	std::valarray<double> _v_g; // Скрость газа
	std::valarray<double> _v_l; // Скрость жидкости
	std::valarray<double> _alpha_g; // Объёмная доля газа
	std::valarray<double> _p; // Давление дисперсной среды


	void PrintCourantNumber(const std::valarray<double>& v_m_intermediate);
	void InitializeGeometryParameters(); // Инициализация параметров скважины
	int CalculateNumberOfPoints();
	void CalculateApproximateMixtureVelocity(std::valarray<double> & v_m_intermediate);
	std::valarray<double> CalculateMixtureVelocityCorrection(const std::valarray<double> & p_corr);
	std::valarray<double> CalculateGasVelocity(const std::valarray<double>& p, const std::valarray<double>& v_m, const std::valarray<double> & alpha_g_previous_iteration);
	std::valarray<double> CalculateLiquidVelocity(const std::valarray<double>& v_m, const std::valarray<double>& alpha_g, const std::valarray<double>& v_g);
	std::valarray<double> CalculateGasVolumeFraction(const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g);
	std::valarray<double> CalculatePressureCorrection(const std::valarray<double> & v_m_intermediate, const std::valarray<double>& p_intermediate);
	double CalculateGasImbalance(const std::valarray<double>& alpha_g_intermediate, const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g_intermediate, const std::valarray<double>& v_l_intermediate);
	void TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d);
	void SimpleAlgorithm();
};

